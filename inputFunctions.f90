MODULE inputFunctions
   
USE constantsMod 
USE inputParametersMod
USE inputFileMod

IMPLICIT NONE 

INTEGER :: k,sgnbdos,j
REAL(DP) :: arg,dt,dz,PHt,expi,expr,args,tt
INTEGER, ALLOCATABLE :: signk(:)
DOUBLE COMPLEX, ALLOCATABLE :: Ht(:)
REAL(DP), ALLOCATABLE :: t(:)
REAL(DP), ALLOCATABLE :: u0t(:,:)
REAL(DP), ALLOCATABLE :: zc(:)
REAL(DP), ALLOCATABLE :: ttau(:)
REAL(DP), ALLOCATABLE :: nChirp(:)
REAL(DP), ALLOCATABLE :: sigma(:)


CONTAINS

!########################
!### stepFunction
!########################

SUBROUTINE stepFunction()			

 dt = (tmax-tmin)/nSteps			!diferencial de dt=0.02
 dz=(zmax-zmin)/nzSteps 			

 allocate(signk(0:nSteps-1))			
 allocate(Ht(0:nSteps-1))			

   do k=0,nSteps-1
     t(k)= tmin + k*dt				!la función t(k) tendra un incremento de 0.02 en cada iteracion
     if(t(k).lt.0) signk(k)=-1
    if(0.le.t(k)) signk(k)= 1
     if(abs(t(k))<=0.5d0) then			
        u0t(1,k)= 1.d0				
     else					
        u0t(1,k)= 0.d0				
     end if					
!    if(0.le.t(k)) then	!signk(k)= 1		
!     arg= -2.d0*pi*freq*t(k)			
!     u0t(1,k)= dcos(arg)			
!     u0t(2,k)= dsin(arg)			
     Ht(k)=cmplx(u0t(1,k),u0t(2,k))
    PHt=dreal(Ht(k))**2+dimag(Ht(k))**2		
!      write(6,*)t(k),Ht(k)
!     write(6,'(3f16.6)')t(k),Ht

   end do ! k

deallocate(signk,Ht)				

END SUBROUTINE stepFunction

! ####################################
! ChirpedgaussianFunction()
!     
! f(t)=exp(-(1+I*C)*t^2/(2t_0))=exp(-(1+I*C)*tau^2/2)
! f(t)=exp(-t^2/(2t_0))*(cos(C*t^2/(2t_0))-I*sin(C*t^2/(2t_0)))
! with tau=t/t_0
!      tau -> t
!####################################
SUBROUTINE gaussianFunction()

 dt = (tmax-tmin)/nSteps
 dz=(zmax-zmin)/nzSteps 

 allocate(signk(0:nSteps-1))
 allocate(Ht(0:nSteps-1))

   do k=0,nSteps-1      
      t(k)= tmin + k*dt
      if(t(k)<0) signk(k)=-1			!k<2499
      if(0<=t(k)) signk(k)= 1			!k>2500
      arg= -(t(k)**2.d0)/(2.d0)				
      expr=cos(arg*Chirp)
      expi=-1*sin(arg*Chirp)			
      u0t(1,k)= expr*exp(arg)
      u0t(2,k)= expi*exp(arg)
      Ht(k)=cmplx(u0t(1,k),u0t(2,k))
      PHt=dreal(Ht(k))**2+dimag(Ht(k))**2
	
!      if(lambda>lambda_dis) sgnbdoa=-1
!      if(lambda<lambda_dis) sgnbdoa=1
      sgnbdos=-1	
      zc(k)=zmin+k*dz*0.01						
      ttau(k)=((1+sgnbdos*Chirp*zc(k)/(T0**2))**2 + (sgnbdos*zc(k)/(T0**2))**2)**0.5      	
      nChirp(k)=Chirp+(1+Chirp**2)*(sgnbdos*zc(k)/(T0**2))
 
!      write(6,*)t(k),Ht(k),PHt
!      write(6,'(4f16.6)')t(k),Ht(k),PHt
   end do ! k

deallocate(signk,Ht)

END SUBROUTINE gaussianFunction

!###############################
!### supergaussianFunction()
!     
! f(t)=exp(-(1+I*C)*(t/(2t_0))^(2*m))=exp(-(1+I*C)*tau^(2*m)/2)
! with tau=t/t_0
!      tau -> t
!###############################
SUBROUTINE supergaussianFunction()

 dt = (tmax-tmin)/nSteps
 dz=(zmax-zmin)/nzSteps 

 allocate(signk(0:nSteps-1))
 allocate(Ht(0:nSteps-1))

   do k=0,nSteps-1      
      t(k)= tmin + k*dt
      if(t(k)<0) signk(k)=-1			!k<2499
      if(0<=t(k)) signk(k)= 1			!k>2500
      arg= -(t(k)**((2.d0)*m))/(2.d0)				
      expr=cos(arg*Chirp)
      expi=-1*sin(arg*Chirp)
      u0t(1,k)= expr*exp(arg)
      u0t(2,k)= expi*exp(arg)
      Ht(k)=cmplx(u0t(1,k),u0t(2,k))
      PHt=dreal(Ht(k))**2+dimag(Ht(k))**2
	
!      write(6,*)t(k),Ht(k),PHt
!      write(6,'(4f16.6)')t(k),Ht(k),PHt
   end do ! k

deallocate(signk,Ht)

END SUBROUTINE supergaussianFunction

!####################################
!### HiperbolicSecantFunction()
!     
! f(t)=2*exp(-I*C*tau^2/2)/(exp(tau)+exp(-tau))
! with tau=t/t_0
!      tau -> t
!####################################
SUBROUTINE HiperSecFunction()

 dt = (tmax-tmin)/nSteps
 dz=(zmax-zmin)/nzSteps 

 allocate(signk(0:nSteps-1))
 allocate(Ht(0:nSteps-1))

   do k=0,nSteps-1      
      t(k)= tmin + k*dt
      if(t(k)<0) signk(k)=-1			!k<2499
      if(0<=t(k)) signk(k)= 1			!k>2500
      arg= -(t(k)**2.d0)/(2.d0)	
      args=1*(t(k))			
      expr=cos(arg*Chirp)
      expi=-sin(arg*Chirp)
      u0t(1,k)=2*expr*(exp(args)+exp(-args))**(-1)
      u0t(2,k)=2*expi*(exp(args)+exp(-args))**(-1)
      Ht(k)=cmplx(u0t(1,k),u0t(2,k))
      PHt=dreal(Ht(k))**2+dimag(Ht(k))**2

!      write(6,*)t(k),Ht(k),PHt
!      write(6,'(4f16.6)')t(k),Ht(k),PHt
   end do ! k

deallocate(signk,Ht)

END SUBROUTINE HiperSecFunction

!####################################

END MODULE inputFunctions
