PROGRAM simple4

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

      INCLUDE 'mpif.h'
      INTEGER status(MPI_STATUS_SIZE)


      CALL mpi_init(ierr)

      CALL mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
      CALL mpi_comm_size(MPI_COMM_WORLD,size,ierr)

      write(greeting,100) my_rank, size

      IF(my_rank.eq.0) THEN
         write(6,*) greeting
         DO partner=1,size-1
         CALL mpi_recv(greeting, 50, MPI_CHARACTER, partner, 1, MPI_COMM_WORLD, status, ierr)
            write(6,*) greeting
         END DO
      else
         CALL solve(my_rank)
         CALL mpi_send(greeting, 50, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, ierr)
      END IF

      IF(my_rank.eq.0) THEN
         write(6,*) 'That is all for now!'
      END IF

      CALL mpi_finalize(ierr)

100   format('Hello World: processor ', I2, ' of ', I2)

END PROGRAM