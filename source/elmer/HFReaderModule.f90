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
!-------------------------------------------------------------------------------
MODULE HFReaderModule
    USE Types
    USE DefUtils
    IMPLICIT NONE
    PUBLIC :: heatFluxOnNodes

    !Interface required to pass function handle through to Elmer
    INTERFACE
        FUNCTION heatFluxOnNodes(Model, n, t) RESULT(hf)
            USE DefUtils
            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: n
            REAL(KIND=dp) :: t, hf, ts
            Logical :: GotIt
            TYPE(ValueList_t), POINTER :: BC
            Character(LEN=255) :: f
            Character(LEN=100) :: nodalPrefix
            Character(LEN=100) :: timeString            
            REAL, ALLOCATABLE :: data(:,:)
            INTEGER :: numLines
            INTEGER :: tsIdx

        END FUNCTION heatFluxOnNodes

        SUBROUTINE ReadCSV(filename, data, numLines)
            CHARACTER(LEN=255) :: filename
            REAL, ALLOCATABLE, INTENT(OUT) :: data(:,:)
            INTEGER, INTENT(OUT) :: numLines
            CHARACTER(LEN=200) :: line
            CHARACTER(LEN=20), DIMENSION(2) :: splitLine
            INTEGER, DIMENSION(2) :: splitPos
            INTEGER :: ioStat, fileUnit, i
            CHARACTER(len=255) :: cwd

        END SUBROUTINE ReadCSV

        SUBROUTINE SplitString(str, splitStr)
            CHARACTER(LEN=*), INTENT(IN) :: str
            CHARACTER(LEN=20), DIMENSION(2), INTENT(OUT) :: splitStr
            INTEGER :: endPos

        END SUBROUTINE SplitString
    END INTERFACE
END MODULE HFReaderModule


! Outside any "scope" the Functions Declared 
! in the Interface above must be implemented
FUNCTION heatFluxOnNodes(Model, n, t) RESULT(hf)
    USE DefUtils
    USE HFReaderModule , except_this_one => heatFluxOnNodes
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    INTEGER :: n
    REAL(KIND=dp) :: t, hf, q_flow
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: BC
    Character(LEN=255) :: f
    Character(LEN=100) :: nodalPrefix
    Character(LEN=100) :: timeString
    REAL, SAVE, ALLOCATABLE :: data(:,:)
    INTEGER :: numLines
    REAL(KIND=dp), SAVE :: tLastRead = -1.0

    !Load the boundary condition with variable for filename
    BC => GetBC()
       
    nodalPrefix = getString(BC, 'nodalHFprefix', GotIt)
    write(timeString, '(F15.9)') t
    f = TRIM(nodalPrefix) // '_' // TRIM(ADJUSTL(timeString)) // '.dat'
    IF(.NOT. GotIt) CALL Fatal('heatFluxOnNodes', 'file: ' //  f // ' not found')
    IF ( .NOT. ASSOCIATED(BC) ) THEN
        CALL FATAL('heatFluxOnNodes','No boundary condition found')
    END IF

    !read file then pass the HF value for each node to g
    !only read the file once for each timestep
    IF(t .NE. tLastRead) THEN
        CALL ReadCSV(f, data, numLines)
        hf = data(n,2)
        tLastRead = t
    ELSE
        hf = data(n,2)
    ENDIF

    IF(hf .LT. 0.0) THEN
        print *, "Found negative HF value..."
        hf = 0.0
    ENDIF
    
    ! Read the heat flux value (q_flow) from the boundary condition
    q_flow = getConstReal(BC, 'q_flow', GotIt)
    IF (.NOT. GotIt) THEN
         print *, "q_flow not found in boundary condition"
         q_flow = 0.0
    END IF    
    
    ! Add constant heat flux to boundary
    hf = hf + q_flow
    
END FUNCTION heatFluxOnNodes

SUBROUTINE ReadCSV(filename, data, numLines)
    USE HFReaderModule, except_this_one => ReadCSV
    CHARACTER(LEN=255) :: filename
    REAL, ALLOCATABLE, INTENT(OUT) :: data(:,:)
    INTEGER, INTENT(OUT) :: numLines
    CHARACTER(LEN=200) :: line
    CHARACTER(LEN=20), DIMENSION(2) :: splitLine
    INTEGER, DIMENSION(2) :: splitPos
    INTEGER :: ioStat, fileUnit, i
    logical :: exists

    print *, "Reading Heat Flux CSV File..."

    ! First pass: Count the number of lines
    fileUnit = 10  ! Arbitrary choice, ensure this unit is not in use elsewhere
    OPEN(UNIT=fileUnit, FILE=filename, ACTION='READ')
    numLines = 0
    DO
        READ(fileUnit, '(A)', IOSTAT=ioStat) line
        IF (ioStat /= 0) EXIT
        numLines = numLines + 1
    END DO
    CLOSE(fileUnit)

    ! Allocate the array
    ALLOCATE(data(numLines, 2))

    ! Second pass: Read the data
    OPEN(UNIT=fileUnit, FILE=filename, ACTION='READ')
    DO i = 1, numLines
        READ(fileUnit, '(A)', IOSTAT=ioStat) line
        IF (ioStat /= 0) EXIT
        CALL SplitString(line, splitLine)
        IF (status /= 0) THEN
            PRINT *, 'Error: Line ', i, ' does not have the correct format.'
            STOP
        END IF        
        READ(splitLine(1), *) data(i, 1)
        READ(splitLine(2), *) data(i, 2)
    END DO
    CLOSE(fileUnit)
END SUBROUTINE ReadCSV

SUBROUTINE SplitString(str, splitStr)
    CHARACTER(LEN=*), INTENT(IN) :: str
    CHARACTER(LEN=20), DIMENSION(2), INTENT(OUT) :: splitStr
    INTEGER :: endPos

    endPos = INDEX(str, ',')
    splitStr(1) = str(1:endPos-1)
    splitStr(2) = str(endPos+1:)

END SUBROUTINE SplitString
