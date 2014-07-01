MODULE inputFunctions
   
USE constantsMod 
USE inputParametersMod
USE inputFileMod

IMPLICIT NONE 

INTEGER :: k
REAL(DP) :: arg,dt,dz,PHt
INTEGER, ALLOCATABLE :: signk(:)
DOUBLE COMPLEX, ALLOCATABLE :: Ht(:)
REAL(DP), ALLOCATABLE :: t(:)
REAL(DP), ALLOCATABLE :: u0t(:,:)

CONTAINS

!########################
!### stepFunction
!########################
SUBROUTINE stepFunction()			!pulso cuadrado

 dt = (tmax-tmin)/nSteps			!diferencial de dt divide (50-(-50))/5000=0.02
 dz=(zmax-zmin)/nzSteps 			!diferencial de dz divide (10-0)/50=0.2

 allocate(signk(0:nSteps-1))			!arreglo sink asignado de dimensión 0 hasta 4999, es decir 5000 elementos
 allocate(Ht(0:nSteps-1))			!arreglo Ht

   do k=0,nSteps-1
     t(k)= tmin + k*dt				!la función t(k) tendra un incremento de 0.02 en cada iteracion
     if(t(k).lt.0) signk(k)=-1
     if(0.le.t(k)) signk(k)= 1
!     if(abs(t(k))<=0.5d0) then			!condicion que se cumple en el intervalo 2475<k<2525
!        u0t(1,k)= 1.d0				!asignamos el valor '1' al vector renglon dentro de este intervalo
!     else					!0<k<2474 U 2536<k<4999
!        u0t(1,k)= 0.d0				!asignamos '0' al vector renglon en todo este intervalo
!     end if					!al final tendremos un vector de 1 renglon y 5000 columnas con '1' 's en 2475<k<2525 y '0' en lo demas
     if(0.le.t(k)) !signk(k)= 1		!para los valores positivos de t(k) es decir para cuando k>2501, el arreglo tendra arreglos de 1's
     arg= -2.d0*pi*freq*t(k)			!argumento de las funciones coseno y seno
     u0t(1,k)= dcos(arg)					!!comentar esta llinea para restaurar cambios
     u0t(2,k)= dsin(arg)			!el segundo renglon de la matriz u0t
     Ht(k)=cmplx(u0t(1,k),u0t(2,k))
    PHt=dreal(Ht(k))**2+dimag(Ht(k))**2				!PHT es el modulo cuadrado 		
!      write(6,*)t(k),Ht(k)
!     write(6,'(3f16.6)')t(k),Ht
   end do ! k

deallocate(signk,Ht)				!Libera la memoria que habia dejado asignada para las funciones y que no fueron ocupadas

END SUBROUTINE stepFunction


!####################################
!### gaussianFunction()
!     
! f(t)=exp(-t^2/(2t_0))=exp(-tau^2/2)
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
      if(t(k)<0) signk(k)=-1			!para los valores de k<2499 se cumple esta condicion y asigna a estos valores de k -1's en el arreglo
      if(0<=t(k)) signk(k)= 1			!para los valores de k>2500 se cumple esta condicion y asigna a estos valores de k 1's en el arreglo
      arg= -(t(k)**2.d0)/(2.d0)					!Doble asterisco significa potencia t(k)²=t(k)**2
      u0t(1,k)= exp(arg)
      u0t(2,k)= 0.0d0
      Ht(k)=cmplx(u0t(1,k),u0t(2,k))
      PHt=dreal(Ht(k))**2+dimag(Ht(k))**2	!PHT es el modulo cuadrado 
!      write(6,*)t(k),Ht(k),PHt
!      write(6,'(4f16.6)')t(k),Ht(k),PHt
   end do ! k

deallocate(signk,Ht)

END SUBROUTINE gaussianFunction

!###############################

END MODULE inputFunctions
