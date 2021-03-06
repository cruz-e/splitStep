!!!######################
MODULE InputParametersMod
!!!######################

  USE constantsMod

  IMPLICIT NONE
  
  INTEGER :: nSteps          ! number steps
  INTEGER :: nzSteps         ! number of z steps
  INTEGER :: m		     ! m for supergaussian pulse

  REAL(DP) :: T0             ! width of the gaussian pulse
  REAL(DP) :: tmin           ! Min time for the FFT
  REAL(DP) :: tmax           ! Max time for the FFT
  REAL(DP) :: Chirp	     ! Chirp grade		
  REAL(DP) :: betaotwo	     ! GVD coefficient	
  REAL(DP) :: betaothree     ! TOD coefficient	
  REAL(DP) :: P0	     ! Laser power
  REAL(DP) :: divsize	     ! Division size preescalator (5 or 10)
  REAL(DP) :: zmin           ! Min value for  z(distance) 
  REAL(DP) :: zmax           ! Max value for z
  REAL(DP) :: freq           ! Freq. for the test function
  REAL(DP) :: tol            ! tolerance
  REAL(DP) :: alpha          ! accounts for fiber losses
  REAL(DP) :: gama           ! 
  
  CHARACTER(LEN=20) :: inputFunction ! input  Function
  CHARACTER(LEN=20) :: case
  CHARACTER(LEN=80) :: paramFile
  
!!!##########################
END MODULE InputParametersMod
!!!##########################
