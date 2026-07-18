!-------------------------------------------------------------------------------
!> File: HFReaderModule.f90
!> Engineer: TL/AR
!> Date: 20240924
!> Reads a comma delimited csv file with columns (nodeId, heatFlux[W/m2]) and 
!> assigns the heat flux to the appropriate node in the Elmer mesh
!> Takes the node number, n, and the timestep, t, and
!> then reads the variable nodalHFprefix from the Boundary Condition
!> in Elmer SIF.  Returns variable, hf, which is the hf for that node
!> 
!> The Boundary Condition in SIF should be structured like this:
!> Boundary Condition 1
!>   Target Boundaries(1) = 1 
!>   Name = "HFreader"
!>   Heat Flux = Variable Time
!>   Real Procedure "HEATLibrary" "heatfluxonnodes"
!>   nodalHFprefix = String nodehf
!> End
!>
!> where nodehf is the prefix for the file that we will read.  this script
!> will append a <timestep>.dat onto the prefix before the read.  So the
!> final filename will look like '<nodalHFprefix>_<timestep>.dat'.  This
!> enables a unique file for the surface flux at each timestep.  The precision
!> of the timestep should be F15.9 (so nanosecond precision).  If using the nodehf
!> prefix the filename at 0.2s would be: 'nodehf_0.200000000.dat'
!> if your timestep is larger than 99999s, you will need to increase the 15
!> in the format specification.
!> 
!> The nodehf.dat file should contain two columns, (nodeId, heat flux [W/m2]),
!> and should be comma delimited.  A few exaple lines:
!>    7.597000000E+03, 2.001488694E+07
!>    7.598000000E+03, 2.086694988E+07
!>    7.599000000E+03, 1.944927053E+07
!>    7.600000000E+03, 0.000000000E+00
!>    7.601000000E+03, 0.000000000E+00
!>    7.602000000E+03, 2.109000231E+07
!>    7.603000000E+03, 0.000000000E+00
!> 
!> Compile like this:  elmerf90 -o HEATLibrary.so HFReaderModule.f90
!> Build grid like this:  ElmerGrid 8 2 <name_of_part> -autoclean -relh 1.0
!> run like this: ElmerSolver case.sif
!
!> note that on some CPU architectures, there may be some floating point
!> discrepancies so you may need to export these env vars before
!> running Elmer:
!> export OPENBLAS_NUM_THREADS=1
!> export OMP_NUM_THREADS=1
!> export OPENBLAS_CORETYPE=Prescott
!-------------------------------------------------------------------------------
MODULE HFReaderCore
    USE Types
    USE DefUtils
    IMPLICIT NONE
    
    !Dynamically allocated arrays to hold the node IDs and heat fluxes from the .csv input file. Split the information into two arrays
    REAL(KIND=dp), ALLOCATABLE :: CSV_Nodes(:), CSV_Fluxes(:)
    !Array to map Elmer's local node index to the index in the CSV arrays
    INTEGER, ALLOCATABLE :: LocalMap(:) 
    
    INTEGER :: numLines = 0                      
    REAL(KIND=dp) :: tLastRead = -1.0_dp         
    LOGICAL :: FirstVisit = .TRUE.               
    LOGICAL :: MapInitialized = .FALSE.          
    CHARACTER(LEN=100) :: nodalPrefix            

CONTAINS

    SUBROUTINE ReadCSV(filename)
        CHARACTER(LEN=*) :: filename
        INTEGER :: ioStat, fileUnit, i
        LOGICAL :: exists
        
        !Check if the file actually exists on the disk before trying to open it
        INQUIRE(FILE=filename, EXIST=exists)
        IF (.NOT. exists) CALL Fatal('ReadCSV', 'File not found: '//TRIM(filename))
        
        OPEN(NEWUNIT=fileUnit, FILE=filename, STATUS='OLD', ACTION='READ')
        
        !Count the number of lines in the file to size our arrays
        numLines = 0
        DO
            READ(fileUnit, *, IOSTAT=ioStat)
            IF (ioStat /= 0) EXIT ! Exit loop when EOF is reached
            numLines = numLines + 1
        END DO
        
        !Clean up old memory and allocate fresh arrays for the new data
        IF (ALLOCATED(CSV_Nodes)) DEALLOCATE(CSV_Nodes, CSV_Fluxes)
        IF (numLines > 0) ALLOCATE(CSV_Nodes(numLines), CSV_Fluxes(numLines))
        
        !Rewind to the start of the file and read the actual data
        REWIND(fileUnit)
        DO i = 1, numLines
            READ(fileUnit, *, IOSTAT=ioStat) CSV_Nodes(i), CSV_Fluxes(i)
        END DO
        
        CLOSE(fileUnit)
    END SUBROUTINE ReadCSV

    SUBROUTINE SortCSV()
        INTEGER :: i, j, gap
        REAL(KIND=dp) :: tmpNode, tmpFlux
        
        gap = numLines / 2
        DO WHILE (gap > 0)
            DO i = gap + 1, numLines
                tmpNode = CSV_Nodes(i)
                tmpFlux = CSV_Fluxes(i)
                j = i
                !Shift elements that are greater than tmpNode to the right
                DO WHILE (j > gap)
                    IF (CSV_Nodes(j-gap) <= tmpNode) EXIT
                    CSV_Nodes(j) = CSV_Nodes(j-gap)
                    CSV_Fluxes(j) = CSV_Fluxes(j-gap)
                    j = j - gap
                END DO
                !Place tmpNode in its correct sorted position
                CSV_Nodes(j) = tmpNode
                CSV_Fluxes(j) = tmpFlux
            END DO
            gap = gap / 2
        END DO
    END SUBROUTINE SortCSV
END MODULE HFReaderCore

FUNCTION heatFluxOnNodes(Model, n, t) RESULT(hf)
    USE DefUtils
    USE HFReaderCore
    IMPLICIT NONE
    
    TYPE(Model_t) :: Model
    INTEGER :: n, i, globalNode, low, high, mid
    REAL(KIND=dp) :: t, hf, q_flow
    LOGICAL :: GotIt
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Element_t), POINTER :: Element
    CHARACTER(LEN=255) :: f
    CHARACTER(LEN=100) :: timeString, prefixNow
    
    ! Resolve BC from the current boundary element so each Target Boundary
    ! gets its own q_flow (-1350 / -3375 / 0), not a value cached from BC #1.
    Element => Model % CurrentElement
    IF (ASSOCIATED(Element)) THEN
        BC => GetBC(Element)
    ELSE
        BC => GetBC()
    END IF
    IF (.NOT. ASSOCIATED(BC)) CALL FATAL('heatFluxOnNodes','No BC found')

    ! Always read per-call: multiple HFreader BCs share this procedure/module.
    prefixNow = GetString(BC, 'nodalHFprefix', GotIt)
    IF (.NOT. GotIt) CALL Fatal('heatFluxOnNodes', 'Keyword nodalHFprefix missing')
    IF (FirstVisit .OR. TRIM(prefixNow) /= TRIM(nodalPrefix)) THEN
        nodalPrefix = prefixNow
        FirstVisit = .FALSE.
        ! Prefix change invalidates CSV cache
        tLastRead = -1.0_dp
        MapInitialized = .FALSE.
    END IF

    q_flow = GetConstReal(BC, 'q_flow', GotIt)
    IF (.NOT. GotIt) q_flow = 0.0_dp
    
    !Only read the file if the simulation time has proceed
    IF (ABS(t - tLastRead) > 1.0e-9_dp) THEN
        WRITE(timeString, '(F15.9)') t
        f = TRIM(nodalPrefix) // '_' // TRIM(ADJUSTL(timeString)) // '.dat'
        
        CALL ReadCSV(f)
        IF (numLines > 1) CALL SortCSV()
        
        tLastRead = t
        MapInitialized = .FALSE. 
    END IF
    
    !This maps the local node index (1 to N) to the row in the sorted CSV array
    IF (.NOT. MapInitialized) THEN
        IF (ALLOCATED(LocalMap)) DEALLOCATE(LocalMap)
        ALLOCATE(LocalMap(Model % Mesh % NumberOfNodes))
        LocalMap = 0
        
        DO i = 1, Model % Mesh % NumberOfNodes
            !Determine the true Global Node ID (important for correct mapping when using MPI parallel runs)
            IF (ASSOCIATED(Model % Mesh % ParallelInfo % GlobalDoFs)) THEN
                globalNode = Model % Mesh % ParallelInfo % GlobalDoFs(i)
            ELSE
                globalNode = i
            END IF
            
            !Find this globalNode in the sorted CSV_Nodes array
            low = 1; high = numLines
            DO WHILE (low <= high)
                mid = (low + high) / 2
                IF (NINT(CSV_Nodes(mid)) == globalNode) THEN
                    LocalMap(i) = mid ! Match found, store the CSV index
                    EXIT
                ELSE IF (NINT(CSV_Nodes(mid)) < globalNode) THEN
                    low = mid + 1
                ELSE
                    high = mid - 1
                END IF
            END DO
        END DO
        MapInitialized = .TRUE.
    END IF
    
    hf = 0.0_dp
    
    !If this node was found in the CSV (LocalMap > 0), grab its flux value
    IF (n >= 1 .AND. n <= Model % Mesh % NumberOfNodes) THEN
        IF (LocalMap(n) > 0) hf = CSV_Fluxes(LocalMap(n))
    END IF
    
    !Clamp negative plasma fluxes to zero (cooling is applied via q_flow below)
    IF (hf < 0.0_dp) hf = 0.0_dp
    
    !Add this BC's constant cooling / background flux from the .sif
    hf = hf + q_flow
    
END FUNCTION heatFluxOnNodes
