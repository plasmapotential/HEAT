!-------------------------------------------------------------------------------
!> File: ReXModule.f90
!> Engineer: TL
!> Date: 20240201
!> 
!> Uses the JMAK equation to calculate the recrystallization fraction given
!> the temperature history.  We adapt the JMAK equation so that the time
!> variable is an effective time, related back to an empirical oven test,
!> rather than time itself.  This enables us to have a time varying temperature
!> history, and still use the Avrami coefficients derived at a reference
!> temperature.  
!>
!> Reads in the Temperature, and calculates the ReX fraction at each mesh node, n
!> requires the Avrami coefficient, k0, the material activation energy, E [eV/atom],
!> the Avrami exponent, n, and the reference temperature at which these coefficients
!> were derived, Tref.  The user should supply these values in the Material section
!> of the SIF, and the ReX calculation will be called once per timestep.  Example
!> SIF:
!>  Material 3
!>  Name = "W"
!>    ...other variables...
!>    ReX = Variable Temperature
!>    Real Procedure "ReXLibrary" "rexonnodes" End
!>    avrami_n = Real 1.0
!>    avrami_k0 = Real 4446601450.0
!>    avrami_E = Real 3.0
!>    avrami_Tref = Real 20.0
!>  End
!>
!>
!> Compile like this:  elmerf90 -o ReXLibrary.so ReXModule.f90
!> Compile like this for a single HEAT shared object:  
!>               elmerf90 -o HEATLibrary.so ReXModule.f90 HFReaderModule.f90
!> Build grid like this:  ElmerGrid 8 2 <name_of_part> -autoclean -relh 1.0
!> run like this: ElmerSolver case.sif
!-------------------------------------------------------------------------------
MODULE ReXCore
    USE Types
    USE DefUtils
    IMPLICIT NONE
    
    !Arrays for reading external initial ReX data
    REAL(KIND=dp), ALLOCATABLE :: CSV_Nodes(:), CSV_ReX(:)
    
    !Global arrays tracking the state of every node in the mesh
    REAL(KIND=dp), ALLOCATABLE :: ReXArray(:)     
    REAL(KIND=dp), ALLOCATABLE :: dataIntegral(:) 
    REAL(KIND=dp), ALLOCATABLE :: dataTprev(:)    
    REAL(KIND=dp), ALLOCATABLE :: dataTimeN(:)    
    REAL(KIND=dp), ALLOCATABLE :: dataAvrami(:)   
    INTEGER :: numLines = 0                       
    
CONTAINS

    SUBROUTINE ReadCSV_ReX(filename)
        CHARACTER(LEN=*) :: filename
        INTEGER :: ioStat, fileUnit, i
        LOGICAL :: exists
        
        INQUIRE(FILE=filename, EXIST=exists)
        IF (.NOT. exists) THEN
            numLines = 0
            RETURN
        END IF
        
        OPEN(NEWUNIT=fileUnit, FILE=filename, STATUS='OLD', ACTION='READ')
        
        !Count lines to allocate appropriately sized arrays
        numLines = 0
        DO
            READ(fileUnit, *, IOSTAT=ioStat)
            IF (ioStat /= 0) EXIT
            numLines = numLines + 1
        END DO
        
        !Allocate memory for the Node IDs and their initial ReX values
        IF (ALLOCATED(CSV_Nodes)) DEALLOCATE(CSV_Nodes, CSV_ReX)
        IF (numLines > 0) ALLOCATE(CSV_Nodes(numLines), CSV_ReX(numLines))
        
        !Rewind and read the data
        REWIND(fileUnit)
        DO i = 1, numLines
            READ(fileUnit, *, IOSTAT=ioStat) CSV_Nodes(i), CSV_ReX(i)
        END DO
        
        CLOSE(fileUnit)
    END SUBROUTINE ReadCSV_ReX

    SUBROUTINE SortCSV_ReX()
        INTEGER :: i, j, gap
        REAL(KIND=dp) :: tmpNode, tmpReX
        
        gap = numLines / 2
        DO WHILE (gap > 0)
            DO i = gap + 1, numLines
                tmpNode = CSV_Nodes(i)
                tmpReX = CSV_ReX(i)
                j = i
                DO WHILE (j > gap)
                    IF (CSV_Nodes(j-gap) <= tmpNode) EXIT
                    CSV_Nodes(j) = CSV_Nodes(j-gap)
                    CSV_ReX(j) = CSV_ReX(j-gap)
                    j = j - gap
                END DO
                CSV_Nodes(j) = tmpNode
                CSV_ReX(j) = tmpReX
            END DO
            gap = gap / 2
        END DO
    END SUBROUTINE SortCSV_ReX
END MODULE ReXCore

FUNCTION ReXOnNodes(Model, n, Temp) RESULT(ReX)
    USE DefUtils
    USE ReXCore
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    INTEGER :: n, i, globalNode, low, high, mid
    REAL(KIND=dp) :: Temp, ReX, c, dt, Time, kB
    REAL(KIND=dp) :: avrami_n, avrami_k0, avrami_E
    TYPE(Variable_t), POINTER :: TimeVar
    TYPE(ValueList_t), POINTER :: mat
    CHARACTER(LEN=255) :: f
    CHARACTER(LEN=100) :: nodalReXPrefix
    LOGICAL :: GotIt
    INTEGER :: TotalNodes
    
    !Boltzmann constant in eV/K (required for the Arrhenius exponential)
    kB = 8.617333e-5_dp 
    TotalNodes = Model % Mesh % NumberOfNodes
    
    !Only read the SIF file for Avrami parameters on the very first pass
    IF (.NOT. ALLOCATED(dataAvrami)) THEN
        mat => GetMaterial()
        IF (.NOT. ASSOCIATED(mat)) CALL FATAL('ReXOnNodes','No material found')
        
        avrami_n  = GetConstReal(mat, 'avrami_n', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'avrami_n missing')
        avrami_k0 = GetConstReal(mat, 'avrami_k0', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'avrami_k0 missing')
        avrami_E  = GetConstReal(mat, 'avrami_E', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'avrami_E missing')
        
        !Store in the module array so we don't have to query the SIF again
        ALLOCATE(dataAvrami(3))
        dataAvrami(1) = avrami_n
        dataAvrami(2) = avrami_k0
        dataAvrami(3) = avrami_E
    ELSE
        !Retrieve cached parameters
        avrami_n  = dataAvrami(1)
        avrami_k0 = dataAvrami(2)
        avrami_E  = dataAvrami(3)
    END IF
    
    !Set up the tracking arrays on the first pass
    IF (.NOT. ALLOCATED(ReXArray)) THEN
        ALLOCATE(ReXArray(TotalNodes))
        ReXArray = 0.0_dp
        
        mat => GetMaterial()
        nodalReXPrefix = GetString(mat, 'nodalReXprefix', GotIt)
        
        !If an initial state file is specified, read and map it to our nodes
        IF (GotIt) THEN
            f = TRIM(nodalReXPrefix) // '.dat'
            CALL ReadCSV_ReX(f)
            
            IF (numLines > 0) THEN
                IF (numLines > 1) CALL SortCSV_ReX()
                
                !Binary search to map global nodes to local mesh nodes
                DO i = 1, TotalNodes
                    IF (ASSOCIATED(Model % Mesh % ParallelInfo % GlobalDoFs)) THEN
                        globalNode = Model % Mesh % ParallelInfo % GlobalDoFs(i)
                    ELSE
                        globalNode = i
                    END IF
                    
                    low = 1; high = numLines
                    DO WHILE (low <= high)
                        mid = (low + high) / 2
                        IF (NINT(CSV_Nodes(mid)) == globalNode) THEN
                            ReXArray(i) = CSV_ReX(mid)
                            EXIT
                        ELSE IF (NINT(CSV_Nodes(mid)) < globalNode) THEN
                            low = mid + 1
                        ELSE
                            high = mid - 1
                        END IF
                    END DO
                END DO
                DEALLOCATE(CSV_Nodes, CSV_ReX)
            ELSE
                CALL Warn('ReXOnNodes', 'Init file ' // TRIM(f) // ' not found. Defaulting ReX to 0.0')
            END IF
        END IF  
        
        !Convert the initial ReX fraction back into an equivalent starting time integral
        ALLOCATE(dataIntegral(TotalNodes))
        dataIntegral = 1.0_dp / avrami_k0 * (-log(1.0_dp - ReXArray))**(1.0_dp / avrami_n)
    END IF
    
    !Initialize Time and Temperature Tracking Arrays ---
    IF (.NOT. ALLOCATED(dataTimeN)) THEN
        ALLOCATE(dataTimeN(TotalNodes))
        dataTimeN = 0.0_dp
    END IF
    IF (.NOT. ALLOCATED(dataTprev)) THEN
        ALLOCATE(dataTprev(TotalNodes))
        dataTprev = Temp
    END IF
    
    !Advance the Integration (Trapezoidal Rule)
    TimeVar => VariableGet(Model % Variables, "Time")
    Time = TimeVar % Values(1)
    
    !It only evaluates on the first pass of a new time step.
    IF (dataTimeN(n) /= Time) THEN
        !Convert timestep size from seconds to hours
        dt = GetTimestepSize() / 3600.0_dp
        
        !Calculate the Arrhenius integrand using the Trapezoidal Rule
        c = dt * ( exp(-avrami_E / (kB * dataTprev(n))) + exp(-avrami_E / (kB * Temp)) ) / 2.0_dp
        
        !Accumulate the integral and update tracking variables
        dataIntegral(n) = dataIntegral(n) + c
        dataTprev(n) = Temp
        dataTimeN(n) = Time
    END IF   
    
    !Apply the standard JMAK equation using the accumulated integral
    ReXArray(n) = 1.0_dp - exp(-(avrami_k0 * dataIntegral(n))**avrami_n )
    
    !Prevent solver crashes if the math yields NaN or Infinity
    IF (ieee_is_nan(ReXArray(n)) .OR. .NOT. ieee_is_finite(ReXArray(n))) THEN
        ReXArray(n) = 0.0_dp
    END IF  
    
    ReX = ReXArray(n)
END FUNCTION ReXOnNodes

SUBROUTINE WriteReX(Model, Solver, dt)
    USE DefUtils
    USE ReXCore
    IMPLICIT NONE
    
    TYPE(Model_t), INTENT(IN) :: Model
    TYPE(Solver_t), INTENT(IN) :: Solver
    REAL(KIND=8), INTENT(IN) :: dt
    LOGICAL :: GotIt
    CHARACTER(LEN=255) :: filename_new
    CHARACTER(LEN=10) :: rankStr
    INTEGER :: i, fu, globalNode
    
    !Do nothing if the ReX array hasn't been initialized yet
    IF (.NOT. ALLOCATED(ReXArray)) RETURN
    
    !Get the output filename from the SIF's section
    filename_new = GetString(Solver % Values, 'Filename', GotIt)
    IF (.NOT. GotIt) CALL Fatal('WriteReX', 'Missing "Filename" in solver section of SIF.')
    
    !MPI Parallel Safety: Append the partition number to the filename so multiple cores don't write over the exact same file.
    IF (ParEnv % PEs > 1) THEN
        WRITE(rankStr, '(I0)') ParEnv % MyPE
        filename_new = TRIM(filename_new) // '.p' // TRIM(rankStr)
    END IF
    
    !Open the file and overwrite any existing data
    OPEN(NEWUNIT=fu, FILE=filename_new, STATUS="REPLACE", ACTION="WRITE")
    
    !Loop through all nodes in this partition and write their global ID and ReX value
    DO i = 1, SIZE(ReXArray)
        IF (ASSOCIATED(Model % Mesh % ParallelInfo % GlobalDoFs)) THEN
            globalNode = Model % Mesh % ParallelInfo % GlobalDoFs(i)
        ELSE
            globalNode = i
        END IF
        WRITE(fu,'(I0,",",E15.9)') globalNode, ReXArray(i)
    END DO
    
    CLOSE(fu)
END SUBROUTINE WriteReX
