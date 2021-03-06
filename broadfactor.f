!###################################################################
!The width of hyperbolic Secant pulses is given by root-mean-square 
!(RMS) width σ defined as:
!
!sigma=[«T^2»-«T»^2]^(1/2) with
!
!«T^n»= (Integral[(T^n * |U(z,T)|^2 * dT);(-infinity,infinity)])/(Integral[(|U(z,T)|^2 * dT);(-infinity,infinity)])
!
!###################################################################

 subroutine rms(tmax,tmin,Np,t,pot,sigos)

      INTEGER :: k
      REAL(DP) :: dt
      REAL(DP), DIMENSION(0:Np-1) :: sigos,pot,t,ttnum1,ttnum2,ttden
	
allocate(sigos(0:Np-1))
allocate(ttnum1(0:Np-1))
allocate(ttnum2(0:Np-1))
allocate(ttden(0:Np-1))

      do k=0,Np-1
	dt=(tmax-tmin)/Np
	ttnum1=ttnum1 + (t(k)**2)*pot(k)*dt
	ttnum2=ttnum2 + t(k)*pot(k)*dt
	ttden=ttden + pot(k)*dt
	sigos=(ttnum1/ttden - (ttnum2/ttden)**2)*0.5

       end do 
deallocate(sigos,ttnum1,ttnum2,ttden)

      end subroutine rms

!###################################################################
!
!Found FWHM and Broadering Factor(tt)
!
!###################################################################

 subroutine fwhmbf(potencia,Np,errpotencia,brfctr)


	INTEGER :: k,lp,lf,rg
	REAL(DP) :: l,lm,vmn,vmx
	REAL(DP), DIMENSION(0,Np-1) :: potencia,errpotencia,epr,epl,brfctr

allocate(errpotencia(0:Np-1))
allocate(epr(0:Np-1))
allocate(epl(0:Np-1))
allocate(brfctr(0:Np-1))

	l=maxval(potencia)						!Location of intensity  max value 
	lp=(dot_product(maxloc(potencia),maxloc(potencia))-1)**0.5	!Location of intensity position max value 
	lm=0.5*l							!Aproximate half-intensity
	errpotencia(:)=abs(potencia(:)-lm)					!Aproximate Location of FWHM if errpot~0
	
	do k=0,Np-1
		if(k<lp) then	
		epl(k)=errpotencia(k)
		epr(k)=0
		end if	
		if(k>lp) then
		epl(k)=0
		epr(k)=errpotencia(k)
		end if
	end do

	vmn=minloc(epl,dim=1,mask=epl>0)			!Look for errpot~0 at the left
	vmx=minloc(epr,dim=1,mask=epr>0)			!Look for errpot~0 at the right
	lf=vmn							!Save positions
	rg=vmx

	if(rg>lf)     brfctr=t(rg)-t(lf)				!Calculate Broadering factor
	if(rg<lf)     brfctr=t(lf)-t(rg)

deallocate(errpotencia,epr,epl,brfctr)

	end subroutine fwhmbf
!##########################################################################
