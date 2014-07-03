c      PROGRAM PROPAGATION
c      

      INTEGER i,j,k,l,m,n,p,nj
      INTEGER iz,jsign,jf,ns,Npf
      PARAMETER (Np=1000)
      DOUBLE PRECISION pi,dn,dae,twopi
      DOUBLE PRECISION sumhk,sumHn
      DOUBLE PRECISION dk,freq,zero,tmin,tmax,t0
      DOUBLE PRECISION z,zstep,zmax,zmin,w,b2,df
      DOUBLE PRECISION fexpr,fexpi,Puzt
      DIMENSION signk(0:Np-1)
      DIMENSION tk(0:Np-1),u0t(2,0:Np-1)
      DIMENSION fn(0:Np-1),dfn(0:Np-1),u0w(2,0:Np-1),u0wf(2,0:Np-1)
      DIMENSION uzt(2,0:Np-1),uzw(2,0:Np-1)
      DIMENSION tk0(0:Np-1),u0tp(2,0:Np-1)
      DOUBLE COMPLEX czero,chtk

      EXTERNAL sdft

c234567
      pi=dacos(-1.d0)
      twopi=2.d0*pi
      czero=(0.d0,0.d0)
      zero=0.0d0
c
c      dk=0.01
      jf=1

c
c     Test  Function
      if(jf.eq.0) then
       tmin= -20.d0
       tmax= 20.0d0
       dk = (tmax-tmin)/Np
        freq=1.d0
       do k=0,Np-1
          tk(k)= tmin + k*dk
c         if(tk(k).lt.0) signk(k)=-1
c          if(0.le.tk(k)) signk(k)= 1
          if(abs(tk(k))<=0.5d0) then
            u0t(1,k)= 1.d0
          else
             u0t(1,k)= 0.d0
          end if
          if(0.le.tk(k)) signk(k)= 1
          dae= -2.d0*pi*freq*tk(k)
c         u0t(1,k)= 1.d0 !dcos(dae)
           u0t(2,k)= 0.0d0 !dsin(dae)
          chtk=cmplx(u0t(1,k),u0t(2,k))
          write(66,'(3f16.6)')tk(k),chtk
       end do ! k
      end if
c     
c       f(t)=exp(-t^2/(2t_0))=exp(-tau^2/2)
c        with tau=t/t_0
c        tau -> t
      if(jf.eq.1) then
       tmin= -20.d0
       tmax= 20.d0
       dk = (tmax-tmin)/Np
        t0=1.d0
       do k=0,Np-1
          tk(k)= tmin + k*dk
          if(tk(k).lt.0) signk(k)=-1
          if(0.le.tk(k)) signk(k)= 1
          dae= -(tk(k)**2.d0)/(2.d0)
          u0t(1,k)= exp(dae)
          u0t(2,k)= 0.0d0
          chtk=cmplx(u0t(1,k),u0t(2,k))
          write(66,'(4f16.6)')tk(k),chtk,(dreal(chtk)**2+dimag(chtk)**2)
       end do ! k
      end if
c
c
         
        do n=0,Np-1
          fn(n)= n/(Np*dk) 
        end do
c
c     Perform Discrete Fourier Transform
      call sdft(1,Np,dk,u0t,u0w)
c
c       do n=0,Np-1 
c         write(67,'(3f16.6)')fn(n),u0w(1,n),u0w(2,n)
c       end do
c      
c     Arrange data in negative and 
c     positive frecuencies in order
c     to plot in the interval (-f_c,fc)
c      H_n=H_{N-n}
        
c        dfn(Np/2)=-fn(Np/2)
        write(67,'(3f16.6)')-fn(Np/2),u0w(1,Np/2),u0w(2,Np/2)
        do n=(Np/2)+1,Np-1
          dfn(n)=-fn(Np/2)+(fn(n)-fn(Np/2))
          write(67,'(3f16.6)')dfn(n),u0w(1,n),u0w(2,n)
        end do
        do n=0,Np/2
          dfn(n)=fn(n)
          write(67,'(3f16.6)')dfn(n),u0w(1,n),u0w(2,n)
        end do
c
c     Perform Inverse Discrete Fourier Transform
c        F^{-1}(H(f)) = h(t)
      call sdft(-1,Np,dk,u0w,u0tp)
c
       do k=0,Np-1 
         write(68,'(3f16.6)')tk(k),u0tp(1,k),u0tp(2,k) 
      end do
c
c     Parseval Theorem
c    sum_{k=0,N-1} |h_k|^2 = (1/N) sum_{n=0,N-1} |H_n|^2
      sumhk=zero
      sumHn=zero
      do j=0,Np-1
         sumhk = sumhk + (u0t(1,j)**2 + u0t(2,j)**2)
         sumHn = sumHn + (u0w(1,j)**2 + u0w(2,j)**2)
      end do
      write(11,*)'sum h_k = ',sumhk*dk
      write(11,*)'sum H_n = ',sumHn/(Np*dk)
ccccc

c  
c      nj=0
c      do n=0,Np-1
c         if (u0w(1,n) .ge. (0.)) then
c            u0wf(1,nj)=u0w(1,n)
c            u0wf(2,nj)=u0w(2,n)
c            write(67,'(3f16.6)')fn(n),u0wf(1,nj),u0wf(2,nj)
c          end if
c         nj=nj+1
c      end do
c      Npf=nj

c  u(zeta,w)= u(0,w)exp{I b_2 w^2 z/2}
c  u(zeta,w)= u(0,w)exp{I To^2 w^2 zeta/2}
c    with zeta=z/L_D  ; To=1
c     zeta-> z
      df=(fn(1)-fn(0))
      zmin=0.d0
      zmax=5.d0
      ns=50
      zstep=(zmax-zmin)/ns
      do iz=1,5   !ns
         z= zmin + iz*zstep*10
         do n=0,Np-1 
            w=twopi*(dfn(n))
c           if(n<=(Np/2)) w=twopi*(dfn(n))
c           if((Np/2)<n) w=twopi*(-dfn(n))
            T0=1.d0
c            u0w(1,n)= sqrt(twopi)*T0*exp(-(T0**2)*(w**2)/2.d0)
c            u0w(2,n)= 0.d0 
          fa= (T0**2)*(w**2)*z/2.d0
          fexpr= cos(fa)
          fexpi= sin(fa)
          uzw(1,n)= u0w(1,n)*fexpr  - u0w(2,n)*fexpi
          uzw(2,n)= u0w(1,n)*fexpi  + u0w(2,n)*fexpr
         end do
         call sdft(-1,Np,dk,uzw,uzt)
       do k=0,Np-1
         Puzt= (uzt(1,k)**2 + uzt(2,k)**2)
         write(12,'(5f16.6)')z,tk(k),uzt(1,k),uzt(2,k),Puzt
       end do

       write(12,*)""
       write(12,*)""

      end do
ccccc
      END !PROPAGATION
ccccccccccccccccccccccccccccc



