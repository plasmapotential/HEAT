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
MODULE ReXModule
    USE Types
    USE DefUtils
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    PUBLIC :: ReXOnNodes

    !Interface required to pass function handle through to Elmer
    INTERFACE
        FUNCTION ReXOnNodes(Model, n, Temp) RESULT(ReX)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(Solver_t) :: Solver
            TYPE(Variable_t), POINTER :: TimeVar
            REAL(KIND=dp) :: Temp, ReX, kB, test
            REAL(KIND=dp) :: avrami_n, avrami_k0, avrami_E, avrami_Tref
            REAL(KIND=dp) :: c1, c2, c3, dt, Tref, Time
            TYPE(ValueList_t), POINTER :: mat
            Logical :: GotIt
            REAL(KIND=dp), ALLOCATABLE :: dataIntegral(:)
            REAL(KIND=dp), ALLOCATABLE :: dataTprev(:)
            REAL(KIND=dp), ALLOCATABLE :: dataTimeN(:)
            REAL(KIND=dp), ALLOCATABLE :: dataAvrami(:)
            INTEGER :: TotalNodes, n
        END FUNCTION ReXOnNodes
    END INTERFACE
END MODULE ReXModule

! Outside any "scope" the Functions Declared 
! in the Interface above must be implemented
FUNCTION ReXOnNodes(Model, n, Temp) RESULT(ReX)
    USE DefUtils
    USE ReXModule , except_this_one => ReXOnNodes
    USE, INTRINSIC :: ieee_arithmetic
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    TYPE(Solver_t) :: Solver
    TYPE(Variable_t), POINTER :: TimeVar
    REAL(KIND=dp) :: Temp, ReX, kB, test
    REAL(KIND=dp) :: c1, c2, c3, dt, Tref, Time
    REAL(KIND=dp) :: avrami_n, avrami_k0, avrami_E, avrami_Tref
    TYPE(ValueList_t), POINTER :: mat
    Logical :: GotIt
    REAL(KIND=dp), SAVE, ALLOCATABLE :: dataIntegral(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE :: dataTprev(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE :: dataTimeN(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE :: dataAvrami(:)
    INTEGER :: TotalNodes, n

    !physics constants
    kB = 8.617333e-5 ! [eV/K]


    !Load the material and avrami values from SIF
    IF (.NOT. ALLOCATED(dataAvrami)) THEN
        mat => GetMaterial()
        !get avrami exponent
        avrami_n = getConstReal(mat, 'avrami_n', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'could not read avrami parameters from SIF')
        IF ( .NOT. ASSOCIATED(mat) ) THEN
            CALL FATAL('ReXNodes','No ReX Solver found')
        END IF
        !get avrami coefficient
        avrami_k0 = getConstReal(mat, 'avrami_k0', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'could not read avrami parameters from SIF')
        IF ( .NOT. ASSOCIATED(mat) ) THEN
            CALL FATAL('ReXNodes','No ReX Solver found')
        END IF
        !get avrami activation energy [eV/atom]
        avrami_E = getConstReal(mat, 'avrami_E', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'could not read avrami parameters from SIF')
        IF ( .NOT. ASSOCIATED(mat) ) THEN
            CALL FATAL('ReXNodes','No ReX Solver found')
        END IF
        !get oven test reference temperature [degC]
        avrami_Tref = getConstReal(mat, 'avrami_Tref', GotIt)
        IF(.NOT. GotIt) CALL Fatal('ReXOnNodes', 'could not read avrami parameters from SIF')
        IF ( .NOT. ASSOCIATED(mat) ) THEN
            CALL FATAL('ReXNodes','No ReX Solver found')
        END IF

        !JMAK equation coefficients
        Tref = avrami_Tref + 273.15
        c1 = exp(-avrami_E / (kB*Tref))
        c2 = exp(avrami_E / (kB*Tref))       

        ALLOCATE(dataAvrami(6))
        dataAvrami(1) = avrami_n
        dataAvrami(2) = avrami_k0
        dataAvrami(3) = avrami_E
        dataAvrami(4) = Tref
        dataAvrami(5) = c1
        dataAvrami(6) = c2

    ELSE
        avrami_n = dataAvrami(1)
        avrami_k0 = dataAvrami(2)
        avrami_E = dataAvrami(3)
        Tref = dataAvrami(4)
        c1 = dataAvrami(5)
        c2 = dataAvrami(6)
        
    END IF

    ! Allocate the ReX array
    IF (.NOT. ALLOCATED(dataIntegral)) THEN
        TotalNodes = Model % Mesh % NumberOfNodes
        print *, "Number of mesh nodes allocated for ReX:", TotalNodes
        ALLOCATE(dataIntegral(TotalNodes))
        dataIntegral = 0.0
        ReX = 0.0
    END IF

    !Allocate the time tracker
    IF (.NOT. ALLOCATED(dataTimeN)) THEN
        TotalNodes = Model % Mesh % NumberOfNodes
        print *, "Number of mesh nodes allocated for TimeN:", TotalNodes
        ALLOCATE(dataTimeN(TotalNodes))
        dataTimeN = 0.0
    END IF

    ! Allocate the Tprev array
    IF (.NOT. ALLOCATED(dataTprev)) THEN
        TotalNodes = Model % Mesh % NumberOfNodes
        print *, "Number of mesh nodes allocated for Tprev:", TotalNodes
        ALLOCATE(dataTprev(TotalNodes))
        dataTprev = 20.0 + 273.15 ! [K]
    END IF

    !get current timestep
    TimeVar => VariableGet( Model % Variables, "Time" )
    Time = TimeVar % Values(1)

    IF (dataTimeN(n) .NE. TIME) THEN
        !get dt in hours
        dt = GetTimestepSize() / 3600.0
        !temperature history integral (using trapz integration for this timestep)
        c3 = dt * ( exp(-avrami_E / (kB * dataTprev(n))) + exp(-avrami_E / (kB * (Temp+273.15))) ) / 2.0
        !update integral
        dataIntegral(n) = dataIntegral(n) + c3

        !For testing, you can print data for a specific node id here
        !IF (n .EQ. 33936) THEN
        !    print *, "True..."
        !    print *, ReX
        !    print *, dataIntegral(n)
        !    print *, (1 - exp(-avrami_k0*c1 * (c2 * dataIntegral(n))**avrami_n ))
        !    print *, dataTprev(n)
        !    print *, Temp + 273.15
        !END IF

        !save Temperature data for next step
        dataTprev(n) = Temp + 273.15
        !Save timestep data for next timestep
        dataTimeN(n) = Time

    END IF   


    !JMAK equation using our effective time
    ReX = (1 - exp(-avrami_k0*c1 * (c2 * dataIntegral(n))**avrami_n ))
    !These values can sometimes be huge or nan.  check it
    IF (ieee_is_nan(ReX) .OR. .NOT. ieee_is_finite(ReX)) THEN
        print *, "NaN or Inf detected in ReX calc...assigning 0..."
        ReX = 0.0
    END IF    


END FUNCTION ReXOnNodes