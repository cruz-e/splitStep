MODULE ConstantsMod

  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = KIND(1.0)            ! single precision
  INTEGER, PARAMETER :: DP = KIND(1.0D0)          ! double precision
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))     ! single precision complex
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0)) ! double precision complex
  
  REAL(DP), PARAMETER :: pi = 3.1415926535897932384626433832795d0

END MODULE ConstantsMod

