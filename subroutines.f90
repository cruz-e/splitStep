MODULE subroutines
   
USE constantsMod 
USE inputParametersMod
USE inputFileMod

IMPLICIT NONE 

CONTAINS

! ###################################
! SUBROUTINE parsevalTheorem(u0t,u0w)
! ###################################
SUBROUTINE parsevalTheorem(u0t,u0w)
   
!IMPLICIT NONE

INTEGER :: k
REAL(DP) :: sumhk,sumHn,tolPT,dt
REAL(DP), DIMENSION(2,0:nSteps-1) :: u0t,u0w

! #########################################
! check of PARSEVAL THEOREM
!
!    sum_{k=0,N-1} |h_k|^2 = (1/N) sum_{n=0,N-1} |H_n|^2
!
! #########################################

  write(6,*)"Checking Parseval Theorem"
  sumhk=0.d0
  sumHn=0.d0
  tolPT=0.000001
  do k=0,nSteps-1
     sumhk = sumhk + (u0t(1,k)**2 + u0t(2,k)**2)
     sumHn = sumHn + (u0w(1,k)**2 + u0w(2,k)**2)
  end do
     dt = (tmax-tmin)/nSteps
     sumhk=sumhk*dt
     sumHn=sumHn/(nSteps*dt)
  if((abs(sumhk-sumHn)<=tolPT)) then
   write(6,*)'sum h_k = ',sumhk
   write(6,*)'sum H_n = ',sumHn
   write(6,*)"PARSEVAL THEOREM IS SATISFIED."
  else
   write(6,*)
   write(6,*)"PARSEVAL THEOREM IS NOT SATISFIED."
   write(6,*)"STOPPING: check code"
   write(6,*)
   STOP
  end if

!#################################
END SUBROUTINE parsevalTheorem
!#################################

! ###################################
! SUBROUTINE arrangeFreqData
!      
!     Arrange data in negative and 
!     positive frecuencies in order
!     to plot in the interval (-f_c,fc)
!     H_n=H_{N-n}
!
! ###################################
SUBROUTINE arrangeFreqData(dt,deltaf)
   
!IMPLICIT NONE

INTEGER :: n
REAL(DP) :: dt
REAL(DP), DIMENSION(0:nSteps-1) :: deltaf
REAL(DP), ALLOCATABLE :: f(:)
! #####################

  allocate(f(0:nSteps-1))
  do n=0,nSteps-1
   f(n)= n/(nSteps*dt)
  end do
!
!  deltaf(nSteps/2)=-f(nSteps/2)
!  write(6,'(3f20.6)')-f(nSteps/2)
!
  do n=(nSteps/2)+1,nSteps-1
   deltaf(n)=-f(nSteps/2)+(f(n)-f(nSteps/2))
!  write(6,'(3f20.6)')deltaf(n)
  end do
  do n=0,nSteps/2
   deltaf(n)=f(n)
 ! write(6,'(3f20.6)')deltaf(n)
  end do
!
 deallocate(f)

!#################################
END SUBROUTINE arrangeFreqData
!#################################

END MODULE subroutines
