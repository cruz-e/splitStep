PROGRAM splitStepDistribuido

      !USE nrtype; USE nr; USE nrutil, ONLY : assert,swap
USE constantsMod
USE inputParametersMod
USE inputFileMod
USE inputFunctions
USE subroutines
USE SplitStep

IMPLICIT NONE
      
      INTEGER ierr,my_rank,size,partner
      CHARACTER*50 greeting
      LOGICAL :: file_exists

      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE)

!************************************************************
! Inicialización MPI
!************************************************************
      CALL mpi_init(ierr)

      CALL mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
      CALL mpi_comm_size(MPI_COMM_WORLD,size,ierr)

!************************************************************
! Nodo maestro
!************************************************************
      IF(my_rank.eq.0) THEN
         write(6,*) '***Corriendo splitStep***'

!************************************************************
! Nodos esclavos. Cada uno resuelve una configuración de splitStep
!************************************************************         
      ELSE
         write (paramFile, "(A19,I1)") "inputParametersFile", my_rank
         INQUIRE(FILE=paramFile, EXIST=file_exists) 

         ! Si no existe el archivo, no se procesa
         IF (.NOT.file_exists) THEN
            write (*, *) "No existe el archivo de entrada ", paramFile
         ELSE
            ! Llamada al procedimiento de resolución del módulo de splitStep
            CALL solve(my_rank)
            !CALL mpi_send(greeting, 50, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
         END IF
      END IF

!************************************************************
! Finalización MPI
!************************************************************
      CALL mpi_finalize(ierr)


END PROGRAM