MODULE inputFileMod
  
  USE ConstantsMod, ONLY : DP, pi
  
  USE DebugMod, ONLY : debug
  
  USE InputParametersMod
  
  LOGICAL :: nSteps_set
  LOGICAL :: nzSteps_set
  LOGICAL :: T0_set
  LOGICAL :: tMin_set 
  LOGICAL :: tMax_set
  LOGICAL :: zmin_set
  LOGICAL :: zmax_set
  LOGICAL :: freq_set
  LOGICAL :: Chirp_set
  LOGICAL :: m_set
  LOGICAL :: divsize_set
  LOGICAL :: tol_set
  LOGICAL :: inputFunction_set
  LOGICAL :: case_set
  LOGICAL :: betaotwo_set
  LOGICAL :: betaothree_set
  LOGICAL :: P0_set
  LOGICAL :: alpha_set 
  LOGICAL :: gama_set
  
CONTAINS


! #########################  
  SUBROUTINE parseParamFile(rank)
    IMPLICIT NONE
    INTEGER, INTENT(IN), OPTIONAL :: rank
    CHARACTER(LEN=255) thisLine
    CHARACTER(LEN=255) string1, string2, tmpString
    CHARACTER(LEN=1), PARAMETER :: commentSignal = '#'
    INTEGER :: ios, lineCounter, lineLength
    INTEGER :: iTmp
    
    IF (debug) WRITE(*,*) "Beginning to parse parameter file."

    IF (present(rank)) THEN
      write (paramFile, "(A19,I1)") "inputParametersFile", rank
    ELSE
      paramFile="inputParametersFile"
    END IF


    OPEN(UNIT=1, FILE=paramFile, FORM='FORMATTED', STATUS='OLD', IOSTAT=ios)	 											
    IF (ios.NE.0) THEN								
       WRITE(6,*) "Error opening file: ", TRIM(paramFile)
       WRITE(6,*) "This should be the parameter file and be the first argument."
       WRITE(6,*) "Stopping"
       STOP
    ELSE
       WRITE(6,'(2A)') " Opening file: ", TRIM(paramFile)			
    END IF
    
    lineCounter = 0
    
    DO
       lineCounter = lineCounter + 1
       READ(UNIT=1,FMT='(A)',IOSTAT=ios) thisLine
       IF (ios.EQ.-1) THEN							
          WRITE(*,*) "End of file found"		
          EXIT
       END IF
       
       IF (debug) WRITE(*,'(A)') TRIM(thisLine)					
       
       thisLine = ADJUSTL(thisLine) !move any leading white space to the back.   
       
       CALL remove_comments(thisLine, commentSignal)
       
       IF (debug) WRITE(*,'(A)') TRIM(thisLine)
       
       lineLength = LEN_TRIM(thisLine)						
       
       IF (lineLength.EQ.0) THEN								
          IF (debug) WRITE(*,'(A,I3)') "Blank line found on line number", lineCounter		
          CYCLE											
       END IF
       
       CALL reduce_white_space(thisLine, LEN(thisLine))
       
       CALL check_for_two_strings(thisLine,lineCounter)
       
       READ(thisLine,FMT=*) string1, string2
       
       IF (TRIM(string1) .EQ. "inputFunction") THEN		
          READ(string2,FMT=*) inputFunction		
          inputFunction_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "case") THEN
          READ(string2,FMT=*) case
          case_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "nSteps") THEN
          READ(string2,FMT='(I6)') nSteps
          nSteps_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "nzSteps") THEN
          READ(string2,FMT='(I6)') nzSteps
          nzSteps_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "T0") THEN
          READ(string2,FMT=*) T0
          T0_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "tmin") THEN
          READ(string2,FMT=*) tmin
          tmin_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "tmax") THEN
          READ(string2,FMT=*) tmax
          tmax_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "zmin") THEN
          READ(string2,FMT=*) zmin
          zmin_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "zmax") THEN
          READ(string2,FMT=*) zmax
          zmax_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "tol") THEN
          READ(string2,FMT=*) tol
          zmax_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "divsize") THEN
          READ(string2,FMT=*) divsize
          divsize_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "freq") THEN
          READ(string2,FMT=*) freq
          freq_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "P0") THEN
          READ(string2,FMT=*) P0
          P0_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "alpha") THEN
	  READ(string2,FMT=*) alpha
          alpha_set = .TRUE. 	
       ELSEIF (TRIM(string1) .EQ. "gama") THEN
	  READ(string2,FMT=*) gama
          gama_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "betaotwo") THEN
          READ(string2,FMT=*) betaotwo
          betaotwo_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "betaothree") THEN
          READ(string2,FMT=*) betaothree
          betaothree_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "Chirp") THEN
          READ(string2,FMT=*) Chirp
          Chirp_set = .TRUE.
       ELSEIF (TRIM(string1) .EQ. "m") THEN
          READ(string2,FMT=*) m
          m_set = .TRUE.
       ELSE
          WRITE(*,*) "Error with keyword ", TRIM(string1)
          STOP 'Keyword not identified'
       END IF
    END DO
    
    CLOSE(1)
    
    IF (.NOT.nSteps_set) THEN
       nSteps = 5000
       WRITE(6,*) "Variable nSteps is not set in input. Using default value of",nSteps
    END IF
    
    IF (.NOT.nzSteps_set) THEN
       nzSteps = 50
       WRITE(6,*) "Variable nzSteps is not set in input. Using default value of ",nzSteps
    END IF
    
    IF (.NOT.T0_set) THEN
       T0 = 1.d0
       WRITE(6,*) "Variable T0 is not set in input. Using default value of ",T0
    END IF
    
    IF (.NOT.tmin_set) THEN
       tmin = -50.d0
       WRITE(6,*) "Variable tmin is not set in input file. Using default value of ",tmin
    END IF
    
    IF (.NOT.tmax_set) THEN
       tmax = 50.d0
       WRITE(6,*) "Variable tmax is not set in input file. Using default value of ",tmax
    END IF
    
    IF (.NOT.divsize_set) THEN
       divsize = 0.5d0
       WRITE(6,*) "Variable divsize is not set in input file. Using default value of ",divsize
    END IF
      
    IF (.NOT.zmin_set) THEN
       zmin = 0.d0
       WRITE(6,*) "Variable zmin is not set in input file. Using default value of ",zmin
    END IF
    
    IF (.NOT.zmax_set) THEN
       zmax = 10.d0
       WRITE(6,*) "Variable zmax is not set in input file. Using default value of ", zmax
    END IF
    
    IF (.NOT.freq_set) THEN
       freq = 1.d0
       WRITE(6,*) "Variable freq is not set in input file. Using default value of ", freq
    END IF

    IF (.NOT.Chirp_set) THEN
       Chirp = 0.d0
       WRITE(6,*) "Variable Chirp is not set in input file. Using default value of ", Chirp
    END IF

    IF (.NOT.alpha_set) THEN
       alpha = 0.d0
       WRITE(6,*) "Variable alpha is not set in input file. Using default value of ", alpha
    END IF

    IF (.NOT.gama_set) THEN
       gama = 1.d0
       WRITE(6,*) "Variable gama is not set in input file. Using default value of ", gama
    END IF
        
    IF (.NOT.betaotwo_set) THEN
       betaotwo = 1.d0
       WRITE(6,*) "Variable betaotwo not set in input. Using default value of ", betaotwo
    END IF
   
    IF (.NOT.betaothree_set) THEN
       betaothree = 0.d0
       WRITE(6,*) "Variable betaothree not set in input. Using default value of ",betaothree
    END IF
   
    IF (.NOT.P0_set) THEN
       P0 = 1.d0
       WRITE(6,*) "Variable P0 not set in input. Using default value of ", P0
    END IF

    IF (.NOT.m_set) THEN
       m = 1
       WRITE(6,*) "Variable m is not set in input file. Using default value of ",m
    END IF

    IF (.NOT.inputFunction_set) THEN
       InputFunction = "Gaussian"
       WRITE(6,*) "Variable inputFunction is not set in input file. Using default value of ",inputFunction
    END IF

    IF (.NOT.case_set) THEN
       case = ""
       WRITE(6,*) "Variable case is not set in input file. Using default value of "
    END IF
    
  END SUBROUTINE parseParamFile
!###############################   



!############################### 
  SUBROUTINE remove_comments(inString,commentSignal)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(INOUT) :: inString				
    CHARACTER(LEN=1), INTENT(IN) :: commentSignal			
    INTEGER :: i, commentLoc
    
    ! remove comments
    commentLoc = INDEX(TRIM(inString),commentSignal)			
										
    										
    !WRITE(6,*) "commentLoc: ", commentLoc
    IF (commentLoc.GT.0) THEN							
       inString = inString(1:commentLoc-1)					
       !WRITE(*,'(A,I3)') "Comment found on line number", lineCounter		
    ELSE IF (commentLoc.LT.0) THEN						
       STOP 'Variable commentLoc is negative!  Error in parser.'
    END IF
  END SUBROUTINE remove_comments
!###############################   
  

!############################### 
  SUBROUTINE reduce_white_space(inString, inStringLength)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(INOUT) :: inString				
    INTEGER, INTENT(IN) :: inStringLength				
    CHARACTER(LEN=inStringLength) :: tmpString
    INTEGER :: i, stringLength, lengthCounter
    CHARACTER(LEN=1) :: thisChar, lastChar, space
    
    !WRITE(6,*) "TRIM(inString): ", TRIM(inString)
    
    space = ' '
    thisChar = inString(1:1)							
    IF (thisChar.EQ.space) THEN
       STOP 'Error in reduce_white_space.  Not expecting a leading blank.'
    END IF
    lengthCounter = 1
    tmpString = thisChar
    stringLength = LEN_TRIM(inString)
    IF (stringLength.GT.1) THEN
       DO i=2,stringLength
          lastChar = thisChar
          thisChar = inString(i:i)
          IF (.NOT.((thisChar.EQ.lastChar).AND.(thisChar.EQ.space))) THEN
             lengthCounter = lengthCounter + 1
             tmpString(1:lengthCounter) = tmpString(1:lengthCounter-1) // thisChar
          END IF
       END DO
    END IF
    inString = tmpString
  END SUBROUTINE reduce_white_space
!############################### 

  
!###############################   
  SUBROUTINE check_for_two_strings(inString,lineNumber)				
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: inString
    INTEGER, INTENT(IN) :: lineNumber
    INTEGER :: spaceLoc
    
    ! Check that there are only two strings in the line.  To do this
    ! first locate the first space, and then check for a second one
    spaceLoc = INDEX(TRIM(inString),' ')  ! first blank location
    IF ((spaceLoc.EQ.0).OR.(spaceLoc.EQ.1)) THEN
       WRITE(6,*) "ERROR: spaceLoc found to be ", spaceLoc
       STOP "Error in parser"
    END IF
    spaceLoc = INDEX( TRIM(inString(spaceLoc+1:)),' ') ! second blank location
    IF (spaceLoc .GT. 0) THEN
       WRITE(*,'(A,I3)') "Found too many arguments on line number ", lineNumber
       STOP "Incorrect input file format"
    END IF
  END SUBROUTINE check_for_two_strings
!############################### 
  
END MODULE inputFileMod
