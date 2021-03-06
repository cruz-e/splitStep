MODULE SplitStep

!USE nrtype; USE nr; USE nrutil, ONLY : assert,swap
USE constantsMod
USE inputParametersMod
USE inputFileMod
USE inputFunctions
USE subroutines

!IMPLICIT NONE
INTEGER :: iz,n,v
INTEGER :: lp,rg,lf,vmn,vmx

REAL(DP) :: z,z1,z2,w,zold,beta,sgnbeta,mud,dtm
REAL(DP) :: num1,num2,den,sig
REAL(DP) :: zold2,zoldt2,zold3,zoldt3
REAL(DP) :: ler,mer,l,lm,merg,leiz,fvm,flp!,lf3
REAL(DP) :: fexpr,fexpi,fv1expr,fv1expi,fv2expr,fv2expi
REAL(DP) :: Puzt,Puzt0
REAL(DP) :: arg1l,arg2l,argnl,argz
REAL(DP) :: LsNL,Lseff,phisNL
REAL(DP) :: z2ld,z3ld,zmid

REAL(DP), ALLOCATABLE :: deltaf(:)
REAL(DP), ALLOCATABLE :: pot(:),pot0(:),uzwnl(:)
REAL(DP), ALLOCATABLE :: errpot(:)
REAL(DP), ALLOCATABLE :: u0w(:,:),uzt(:,:),uzw(:,:),v1uzw(:,:),v1uzt(:,:),v2uzt(:,:)
REAL(DP), ALLOCATABLE :: v2uzw(:,:),uz0t(:,:),uz0w(:,:)
REAL(DP), ALLOCATABLE :: u0tp(:,:)
REAL(DP), ALLOCATABLE :: erplf(:)
REAL(DP), ALLOCATABLE :: erprg(:)
REAL(DP), ALLOCATABLE :: ttnum1(:)
REAL(DP), ALLOCATABLE :: ttnum2(:)
REAL(DP), ALLOCATABLE :: ttden(:)
REAL(DP), ALLOCATABLE :: sigos(:)

 CHARACTER(LEN=50)  :: outputPowerFile

 CONTAINS

!###################################################################
! BEGIN
!###################################################################

SUBROUTINE solve(rank) 
 INTEGER, INTENT(IN), OPTIONAL :: rank
 write(6,*)
 write(6,*)"START RUNNING splitStep.f90"
 write(6,*)

 IF (present(rank)) THEN
      CALL parseParamFile(rank)
 ELSE
      CALL parseParamFile
 END IF
 
 write(6,*)'case ', case
 write(6,*)'inputFunction ', inputFunction

  IF (present(rank)) THEN
    write (outputPowerFile, "(A21,I1)") TRIM(case)//'_power-', rank
  ELSE
      outputPowerFile=TRIM(case)//'_power'
 END IF
 OPEN(UNIT=12,FILE=outputPowerFile,FORM='FORMATTED', STATUS='UNKNOWN')

! ###############
! Input Functions
 allocate(u0t(2,0:nSteps-1),t(0:nSteps-1),zc(0:nSteps-1))	
 allocate(u0w(2,0:nSteps-1))
 allocate(ttau(0:nSteps-1))
 allocate(nChirp(0:nSteps-1))
 allocate(sigma(0:nSteps-1))
 allocate(pot(0:nSteps-1))
 allocate(pot0(0:nSteps-1))
 allocate(errpot(0:nSteps-1))
 allocate(erplf(0:nSteps-1))
 allocate(erprg(0:nSteps-1))
allocate(ttnum1(0:nSteps-1))
allocate(sigos(0:nSteps-1))
allocate(ttnum2(0:nSteps-1))
allocate(ttden(0:nSteps-1))
 
 SELECT CASE(inputFunction)
  CASE("stepFunction")
   call stepFunction()
  CASE("Gaussian")
   call gaussianFunction()
  CASE("SuperGaussian")
   call supergaussianFunction()
  CASE("HiperSecFunction")
   call HiperSecFunction()
 END SELECT

!### Perform Discrete Fourier Transform
  write(6,*)"Performing Discrete Fourier Transform"
  call sdft(1,nSteps,dt,u0t,u0w)
! call four1_dp(Ht,1)
!
!### Check Parseval Theorem
 call parsevalTheorem(u0t,u0w)
!###

 allocate(deltaf(0:nSteps-1))
!### Arrange Frequency Data
 call arrangeFreqData(dt,deltaf)
!###

 allocate(u0tp(2,0:nSteps-1))
!### Perform Inverse Discrete Fourier Transform
!        F^{-1}(H(f)) = h(t)
  write(6,*)"Performing Inverse Discrete Fourier Transform"
  call sdft(-1,nSteps,dt,u0w,u0tp)
  do k=0,nSteps-1 
     t(k)= tmin + k*dt
!     write(6,'(3f16.6)')t(k),u0tp(1,k),u0tp(2,k) 
  end do
!###
!### SPLIT-STEP METHOD
!  CASE:  
!  u(zeta,w)= u(0,w)exp{I b_2 w^2 z/2}
!  u(zeta,w)= u(0,w)exp{I To^2 w^2 zeta/2}
!    with zeta=z/L_D  ; To=1
!     zeta-> z
  dz=.1
!write(6,*)'nzSteps ',nzSteps
!write(6,*)'iz, z1, z2'
!############################################ Calculus of U(0,T) is neccesary to get self-phase modulation
   allocate(uz0w(2,0:nSteps-1),uz0t(2,0:nSteps-1))
           do n=0,nSteps-1 
    	      argz=0
              fexpr= cos(argz)
              fexpi= sin(argz)
              uz0w(1,n)= u0w(1,n)*fexpr  - u0w(2,n)*fexpi
              uz0w(2,n)= u0w(1,n)*fexpi  + u0w(2,n)*fexpr
	   end do
         call sdft(-1,nSteps,dt,uz0w,uz0t)
         do k=0,nSteps-1
            Puzt0= (uz0t(1,k)**2 + uz0t(2,k)**2)
	    pot0(k)=Puzt0						!|U(0,T)|
         end do
    deallocate(uz0t,uz0w)
!############################################ A partir de aqui se prueba el cambio (sustituir info archivo sincambios)
   do iz=0,10	
   allocate(uzw(2,0:nSteps-1),uzt(2,0:nSteps-1),v1uzw(2,0:nSteps-1),v2uzw(2,0:nSteps-1),uzwnl(0:nSteps-1))
   allocate(v1uzt(2,0:nSteps-1),v2uzt(2,0:nSteps-1))
	   z= zmin + iz*dz*divsize
	   zmid= dz*divsize/2
	   z2ld=z*abs(betaotwo)/T0**2
	   z3ld=z*abs(betaothree)/T0**3	
	   z1= zmin + iz*dz*divsize
   	   zold2=z1*abs(betaotwo)/T0**2
           zold3=z1*abs(betaothree)/T0**3				
	   z2= zmin + iz*dz*divsize+zmid
   	   zoldt2=z2*abs(betaotwo)/T0**2
           zoldt3=z2*abs(betaothree)/T0**3
!########################################### Non Linear part contribution (NLPC)
!	 do k=0, nSteps-1
!	     LsNL=(gama*P0)**(-1)
!	     if (alpha==0) Lseff=z
!	     if (alpha/=0) Lseff=(1-exp(-alpha*z))/alpha.
!	     uzwnl(k)=exp(-alpha*z)*pot0(k)/LsNL
!	     argnl=argnl+uzwnl(k)*dz
!         end do
!###########################################
!      write(6,*)iz,z1,z2		
	   do n=0,nSteps-1 
              w=2.d0*pi*(deltaf(n))
!             if(n<=(nSteps/2)) w=2.d0*pi*(deltaf(n))
!             if((nSteps/2)<n) w=2.d0*pi*(-deltaf(n))
!             u0w(1,n)= sqrt(2.d0*pi)*T0*exp(-(T0**2)*(w**2)/2.d0)
!             u0w(2,n)= 0.d0 
              arg1l= (T0**2)*(w**2)*zold2/2.d0+(T0**3)*(w**3)*zold3/6.d0		
              arg2l= (T0**2)*(w**2)*zoldt2/2.d0+(T0**3)*(w**3)*zoldt3/6.d0
              argnl=0 									!Integration over z->z+h
	      arg=(T0**2)*(w**2)*z2ld/2.d0+(T0**3)*(w**3)*z3ld/6.d0		
!    	      arg=arg1l+argnl+arg2l	

              fexpr= cos(arg)
              fexpi= sin(arg)
              fv1expr= cos(arg1l)
              fv1expi= sin(arg1l)
              fv2expr= cos(arg2l)
              fv2expi= sin(arg2l)
	      v1uzw(1,n)= u0w(1,n)*fv1expr - u0w(2,n)*fv1expi
	      v1uzw(2,n)= u0w(1,n)*fv1expi + u0w(2,n)*fv1expr
	      v2uzw(1,n)= fv2expr
	      v2uzw(2,n)= fv2expi  
!	      v2uzw(1,n)= u0w(1,n)*fv2expr  - u0w(2,n)*fv2expi
!	      v2uzw(2,n)= u0w(1,n)*fv2expi  + u0w(2,n)*fv2expr	
             uzw(1,n)=u0w(1,n)*fexpr  - u0w(2,n)*fexpi
             uzw(2,n)=u0w(1,n)*fexpi  + u0w(2,n)*fexpr
!             uzw(1,n)=v1uzw(1,n)*v2uzw(1,n) - v1uzw(2,n)*v2uzw(2,n) 			
!             uzw(2,n)=v1uzw(1,n)*v2uzw(2,n) + v1uzw(2,n)*v2uzw(1,n) 			
	   end do
         call sdft(-1,nSteps,dt,uzw,uzt)
!         call sdft(-1,nSteps,dt,v1uzw,v1uzt)
!         call sdft(-1,nSteps,dt,v2uzw,v2uzt)
         do k=0,nSteps-1
!             uzt(1,k)=v1uzt(1,k)*v2uzt(1,k) - v1uzt(2,k)*v2uzt(2,k) 	
!             uzt(2,k)=v1uzt(1,k)*v2uzt(2,k) + v1uzt(2,k)*v2uzt(1,k) 	
!            uzt(1,k)=v1uzt(1,k)*v2uzt(1,k)
!            uzt(2,k)=v1uzt(2,k)*v2uzt(2,k)
            Puzt= (uzt(1,k)**2 + uzt(2,k)**2)
	    pot(k)=Puzt
!            if(iz==0) then
!	      Puzt0=Puzt 
!	      pot0(k)=Puzt0 
!	    end if
         end do
!############################################ A partir de aqui se termina el cambio
!###################################################################
!Found FWHM and Broadering Factor(tt)
!###################################################################
	l=maxval(pot)								!Location of intensity  max value 
	lp=(dot_product(maxloc(pot),maxloc(pot))-1)**0.5			!Location of intensity position max value 
	lm=0.5*l
	errpot(:)=abs(pot(:)-lm)						!Aproximate Location of FWHM if errpot~0
	do k=0,nSteps-1
	   if(k<lp) then	
	      erplf(k)=errpot(k)
	      erprg(k)=0
	   end if	
	   if(k>lp) then
	      erplf(k)=0
	      erprg(k)=errpot(k)
	   end if
	end do
	vmn=minloc(erplf,dim=1,mask=erplf>0)					!Look for errpot~0 at the left
	vmx=minloc(erprg,dim=1,mask=erprg>0)					!Look for errpot~0 at the right
	lf=vmn									!Save positions
	rg=vmx
        tt=t(rg)-t(lf)								!Calculate Broadering factor
!###################################################################
!(RMS) width σ defined as:
!
!sigma=[«T^2»-«T»^2]^(1/2) with
!«T^n»= (Integral[(T^n * |U(z,T)|^2 * dT);(-infinity,infinity)])/(Integral[(|U(z,T)|^2 * dT);(-infinity,infinity)])
!###################################################################
      do k=0,nSteps-1
	dtm=(tmax-tmin)/nSteps
	ttnum1(k)=(t(k)**2)*pot(k)*dtm
	num1=num1+ttnum1(k)
	ttnum2(k)=t(k)*pot(k)*dtm
	num2=num2+ttnum2(k)
	ttden(k)=pot(k)*dtm
	den=den+ttden(k)
	sig=(num1/den-(num2/den)**2)**(0.5)
      end do 
!###################################################################
	do k=0,nSteps-1
!	        write(12,'(7f16.6)')z,zold,t(k),zc(k),ttau(k),nChirp(k),tt		!,errpot(k)
!	        write(12,'(8f16.6)')z,t(k),zold3,pot(k),num1,num2,den,sig!,ttau(k)!,zc(k),ttau(k),nChirp(k)		!,errpot(k)
!	        write(12,'(9f16.6)')z,t(k),v1uzw(1,k),v1uzw(2,k),uzw(1,k),uzw(2,k),uzt(1,k),uzt(2,k),pot(k)
	end do

     write(12,*)""
     write(12,*)""
    deallocate(uzw,uzt,v1uzw,v2uzw,uzwnl,v1uzt,v2uzt)
  end do

  write(6,*)
  write(6,*)'Power output-file: ',outputPowerFile
  write(6,*)

 deallocate(u0t,deltaf,t,zc,ttau,nChirp,sigos,pot,pot0)	
 deallocate(errpot,erplf,erprg)	
 deallocate(ttnum1,ttnum2,ttden)

!###################################################################
 write(6,*)""
 write(6,*)"END PROGRAM splitStep"
!###################################################################
!END PROGRAM splitStep
!###################################################################
END SUBROUTINE solve

END MODULE splitStep