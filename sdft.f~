!!   Subroutine
!  Perform Discrete Fourier Transform
!
      subroutine sdft(jsign,Np,dk,shtk,sHfn)

      INTEGER :: i,j,k,l,m,n,p,jsign
      DOUBLE PRECISION :: pi,arg,dk
      DOUBLE PRECISION, DIMENSION(0:Np-1) :: stk,sfn shtk
      DOUBLE PRECISION, DIMENSION(2,0:Np-1) :: shtk,sHfn
      DOUBLE COMPLEX :: czero,chtk,ce,csH,cft

      
      pi=dacos(-1.d0)
      czero=(0.d0,0.d0)

        do l=0,Np-1 
        csH=czero
        do p=0,Np-1
          if(jsign.eq.1) then
             n=l 
             k=p
          end if
          if(jsign.eq.-1) then
             k=l
             n=p
          end if   
          chtk=cmplx(shtk(1,p),shtk(2,p))
!          if(jsign.eq.-1) write(66,*)chtk
          if(n.eq.1) write(66,*)k,shtk(1,k),shtk(2,k)           
!          arg=jsign*2.d0*pi*sfn(n)*stk(k)
          arg=jsign*2.d0*pi*n*k/Np
          ce=cmplx(cos(arg),sin(arg))
          csH= csH + chtk*ce
!          write(65,*)l,k,dreal(csH),dimag(csH)
!          if(jsign.eq.-1) write(66,*)chtk,ce,csH
        end do
!          write(65,*)""
!          write(65,*)""
          if(jsign.eq.1) cft= dk*csH
          if(jsign.eq.-1) cft= csH/(Np*dk)
          sHfn(1,l)=dreal(cft)
          sHfn(2,l)=dimag(cft)
!         write(67,'(3f16.6)')sfn(n),dreal(cft),dimag(cft)
!          write(67,'(3f16.6)')sfn(n),sHfn(1,n),sHfn(2,n)
       end do 
!
      continue
      end !subroutine
