! This program is the Fortran 90 version of gafortomp.f
!*************************************************************************
MODULE gafort_globals
IMPLICIT NONE
INTEGER, PARAMETER :: DP=KIND(1.0D0)
INTEGER, PARAMETER :: K4=KIND(1)
INTEGER, PARAMETER :: OMP_LOCK_KIND = selected_int_kind(18) 
INTEGER(K4), PARAMETER :: indmax=400000, nchrmax=450, nparmax=30
REAL(DP), DIMENSION(nparmax,indmax) :: parent, child
REAL(DP), DIMENSION(indmax) :: fitness
INTEGER(K4), DIMENSION(nparmax) :: nposibl, nichflg
INTEGER(K4), DIMENSION(nchrmax,indmax) :: iparent, ichild
REAL(DP), DIMENSION(nparmax) :: g0, g1, parmax, parmin, pardel
INTEGER(K4), DIMENSION(nparmax) :: ig2
INTEGER(KIND=OMP_LOCK_KIND), DIMENSION(indmax) :: lck
REAL(DP) :: pcross,pmutate,pcreep
INTEGER(K4) :: maxgen,idum,irestrt,itourny,ielite,icreep,iunifrm, &
               iniche,iskip,iend,nchild,microga,kountmx,npopsiz, &
               nowrite,nparam,nchrome
INTEGER(K4) :: numthreads
INTEGER*8   :: maincnt,evaloutcnt,childcnt,mutatecnt,newgencnt, &
               crosovrcnt,selectcnt,shufflecnt
END MODULE gafort_globals
!*************************************************************************
MODULE ran_state
USE gafort_globals, ONLY : K4,DP
IMPLICIT NONE
INTEGER(K4), PARAMETER :: madim=55
TYPE ranstate
     REAL(DP) :: ma(madim),mj
     INTEGER(K4) :: iff,inext,inextp,iduma
     INTEGER(K4) :: padding(128)
END TYPE ranstate
TYPE(ranstate), ALLOCATABLE :: procstates(:)
INTEGER(K4) :: lenran
!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.                       
REAL(DP), PARAMETER :: mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig
END MODULE ran_state
!*************************************************************************
!
! This is a specially modified version of ga164.f for the SPECOMP2001 
! benchmark suite.  There is some discusion of this GA in the ReadMe
! file.  Note that this version has the restart feature disabled and
! the input file no longer a NAMELIST input.
!
! To download the original version 1.6.4, go to website:
!    <http://cuaerospace.com/carroll/ga.html>
!
! This update on 4/14/99.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! David L. Carroll
! CU Aerospace
! 2004 South Wright Street Extended
! Urbana, IL  61802
!
! e-mail: carroll@cuaerospace.com
! Phone:  217-333-8274
! fax:    217-244-7757
!
! This genetic algorithm (GA) driver is free for public use.  My only
! request is that the user reference and/or acknowledge the use of this
! driver in any papers/reports/articles which have results obtained
! from the use of this driver.  I would also appreciate a copy of such
! papers/articles/reports, or at least an e-mail message with the
! reference so I can get a copy.  Thanks.
!
! This program is a FORTRAN version of a genetic algorithm driver.
! This code initializes a random sample of individuals with different
! parameters to be optimized using the genetic algorithm approach, i.e.
! evolution via survival of the fittest.  The selection scheme used is
! tournament selection with a shuffling technique for choosing random
! pairs for mating.  The routine includes binary coding for the
! individuals, jump mutation, creep mutation, and the option for
! single-point or uniform crossover.  Niching (sharing) and an option
! for the number of children per pair of parents has been added.
! An option to use a micro-GA was recently added.
!
! For companies wishing to link this GA driver with an existing code,
! I am available for some consulting work.  Regardless, I suggest
! altering this code as little as possible to make future updates
! easier to incorporate.
!
! Any users new to the GA world are encouraged to read David Goldberg's
! "Genetic Algorithms in Search, Optimization and Machine Learning,"
! Addison-Wesley, 1989.
!
! Other associated files are:  ga.inp (gafort.in)
!                              ga.out (gafort.out)
!                              ga.restart (not used in SPEC version)
!                              params.f (not used in SPEC version)
!                              ReadMe
!
! I have provided a sample subroutine "func", but ultimately
! the user must supply this subroutine "func" which should be your
! cost function.  You should be able to run the code with the
! sample subroutine "func" and the provided ga.inp file and obtain
! the optimal function value of 1.0 at generation 97 with the
! uniform crossover micro-GA enabled (this is 485 function evaluations).
!
! The code is presently set for a maximum population size of 200,
! 30 chromosomes (binary bits) and 8 parameters.  These values can be
! changed in params.f as appropriate for your problem.  Correspondingly
! you will have to change a few 'write' and 'format' statements if you
! change nchrome and/or nparam.  In particular, if you change nchrome
! and/or nparam, then you should change the 'format' statement numbers
! 1050, 1075, 1275, and 1500 (see ReadMe file).
!
! Please feel free to contact me with questions, comments, or errors
! (hopefully none of latter).
!
! Disclaimer:  this program is not guaranteed to be free of error
! (although it is believed to be free of error), therefore it should
! not be relied on for solving problems where an error could result in
! injury or loss.  If this code is used for such solutions, it is
! entirely at the user's risk and the author disclaims all liability.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input variable definitions:
!
! icreep   = 0 for no creep mutations
!          = 1 for creep mutations; creep mutations are recommended.
! idum     The initial random number seed for the GA run.  Must equal
!          a negative integer, e.g. idum=-1000.
! ielite   = 0 for no elitism (best individual not necessarily
!              replicated from one generation to the next).
!          = 1 for elitism to be invoked (best individual replicated
!              into next generation); elitism is recommended.
! iend         = 0 for normal GA run (this is standard).
!          = number of last population member to be looked at in a set
!            of individuals.  Setting iend-0 is only used for debugging
!            purposes and is commonly used in conjunction with iskip.
! iniche   = 0 for no niching
!          = 1 for niching; niching is recommended.
! irestrt  = 0 for a new GA run, or for a single function evaluation
!          = 1 for a restart continuation of a GA run.
! iskip    = 0 for normal GA run (this is standard).
!          = number in population to look at a specific individual or
!            set of individuals.  Setting iskip-0 is only used for
!            debugging purposes.
! itourny  No longer used.  The GA is presently set up for only
!          tournament selection.
! iunifrm  = 0 for single-point crossover
!          = 1 for uniform crossover; uniform crossover is recommended.
! kountmx  = the maximum value of kount before a new restart file is
!            written; presently set to write every fifth generation.
!            Increasing this value will reduce I/O time requirements
!            and reduce wear and tear on your storage device
! maxgen   The maximum number of generations to run by the GA.
!          For a single function evaluation, set equal to 1.
! microga  = 0 for normal conventional GA operation
!          = 1 for micro-GA operation (this will automatically reset
!            some of the other input flags).  I recommend using
!            npopsiz=5 when microga=1.
! nchild   = 1 for one child per pair of parents (this is what I
!              typically use). This code is written to use 1 child.
!          = 2 for two children per pair of parents (2 is more common
!              in GA work).
! nichflg  = array of 1/0 flags for whether or not niching occurs on
!            a particular parameter.  Set to 0 for no niching on
!            a parameter, set to 1 for niching to operate on parameter.
!            The default value is 1, but the implementation of niching
!            is still controlled by the flag iniche.
! nowrite  = 0 to write detailed mutation and parameter adjustments
!          = 1 to not write detailed mutation and parameter adjustments
! nparam   Number of parameters (groups of bits) of each individual.
!          Make sure that nparam matches the number of values in the
!          parmin, parmax and nposibl input arrays.
! npopsiz  The population size of a GA run (typically 100 works well).
!          For a single calculation, set equal to 1.
! nposibl  = array of integer number of possibilities per parameter.
!            For optimal code efficiency set nposibl=2**n, i.e. 2, 4,
!            8, 16, 32, 64, etc.
! parmax   = array of the maximum allowed values of the parameters
! parmin   = array of the minimum allowed values of the parameters
! pcreep   The creep mutation probability.  Typically set this
!          = (nchrome/nparam)/npopsiz.
! pcross   The crossover probability.  For single-point crossover, a
!          value of 0.6 or 0.7 is recommended.  For uniform crossover,
!          a value of 0.5 is suggested.
! pmutate  The jump mutation probability.  Typically set = 1/npopsiz.
!
!
! For single function evaluations, set npopsiz=1, maxgen=1, & irestrt=0.
!
! My favorite initial choices of GA parameters are:
!    microga=1, npopsiz=5, inunifrm=1, maxgen=100
!    microga=1, npopsiz=5, inunifrm=0, maxgen=100
! I generally get good performance with both the uniform and single-
! point crossover micro-GA.
!
! For those wishing to use the more conventional GA techniques,
! my old favorite choice of GA parameters was:
!    iunifrm=1, iniche=1, ielite=1, itourny=1, nchild=1
! For most problems I have dealt with, I get good performance using
!    npopsiz=100, pcross=0.5, pmutate=0.01, pcreep=0.02, maxgen=26
! or
!    npopsiz= 50, pcross=0.5, pmutate=0.02, pcreep=0.04, maxgen=51
!
! Any negative integer for idum should work.  I typically arbitrarily
! choose idum=-10000 or -20000.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code variable definitions (those not defined above):
!
! best     = the best fitness of the generation
! child    = the floating point parameter array of the children
! cpu      = cpu time of the calculation
! creep    = +1 or -1, indicates which direction parameter creeps
! delta    = del/nparam
! diffrac  = fraction of total number of bits which are different
!            between the best and the rest of the micro-GA population.
!            Population convergence arbitrarily set as diffrac<0.05.
! fbar     = average fitness of population
! fitness  = array of fitnesses of the parents
! fitsum   = sum of the fitnesses of the parents
! g0       = lower bound values of the parameter array to be optimized.
!            The number of parameters in the array should match the
!            dimension set in the above parameter statement.
! g1       = the increment by which the parameter array is increased
!            from the lower bound values in the g0 array.  The minimum
!            parameter value is g0 and the maximum parameter value
!            equals g0+g1*(2**g2-1), i.e. g1 is the incremental value
!            between min and max.
! ig2      = array of the number of bits per parameter, i.e. the number
!            of possible values per parameter.  For example, ig2=2 is
!            equivalent to 4 (=2**2) possibilities, ig2=4 is equivalent
!            to 16 (=2**4) possibilities.
! ig2sum   = sum of the number of possibilities of ig2 array
! ibest    = binary array of chromosomes of the best individual
! ichild   = binary array of chromosomes of the children
! icount   = counter of number of different bits between best
!            individual and other members of micro-GA population
! icross   = the crossover point in single-point crossover
! indmax   = maximum # of individuals allowed, i.e. max population size
! iparent  = binary array of chromosomes of the parents
! istart   = the generation to be started from
! jbest    = the member in the population with the best fitness
! jelite   = a counter which tracks the number of bits of an individual
!             which match those of the best individual
! jend     = used in conjunction with iend for debugging
! jstart   = used in conjunction with iskip for debugging
! kount    = a counter which controls how frequently the restart
!            file is written
! kelite   = kelite set to unity when jelite=nchrome, indicates that
!            the best parent was replicated amongst the children
! mate1    = the number of the population member chosen as mate1
! mate2    = the number of the population member chosen as mate2
! nchrmax  = maximum # of chromosomes (binary bits) per individual
! nchrome  = number of chromosomes (binary bits) of each individual
! ncreep   = # of creep mutations which occurred during reproduction
! nmutate  = # of jump mutations which occurred during reproduction
! nparmax  = maximum # of parameters which the chromosomes make up
! paramav  = the average of each parameter in the population
! paramsm  = the sum of each parameter in the population
! parent   = the floating point parameter array of the parents
! pardel   = array of the difference between parmax and parmin
! rand     = the value of the current random number
! npossum  = sum of the number of possible values of all parameters
! time0    = clock time at start of run
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutines:
! ____________
!
! code     = Codes floating point value to binary string.
! crosovr  = Performs crossover (single-point or uniform). This is inlined.
!	     It is no longer called from main generation loop (compiler issue!).
! decode   = Decodes binary string to floating point value.
! evalout  = Evaluates the fitness of each individual and outputs
!            generational information to the 'ga.out' file.
! func     = The function which is being evaluated.
! gamicro  = Implements the micro-GA technique.
! input    = Inputs information from the 'ga.inp' file.
! initial  = Program initialization and inputs information from the
!            'ga.restart' file.
! mutate   = Performs mutation (jump and/or creep).
! newgen   = Writes child array back into parent array for new
!            generation; also checks to see if best individual was
!            replicated (elitism).
! niche    = Performs niching (sharing) on population.
! possibl  = Checks to see if decoded binary string falls within
!            specified range of parmin and parmax.
! ran3     = The random number generator.
! restart  = Writes the 'ga.restart' file.
! selectn  = Performs selection; tournament selection is the only
!            option in this version of the code. This is inlined. No
!	     longer called from the main generation loop (compiler issue!).
! shuffle  = Shuffles the population randomly for selection.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM gafortran
USE gafort_globals
USE ran_state
IMPLICIT NONE
INTEGER(K4) :: kount, i, istart, ig2sum, npossum, ncross, ipick, ifirst, &
               isecond, iother, mate1, mate2, itemp, n, my_cpu_id, j, ii, &
               icross,k,offset,ithird,iforth
REAL(DP) :: fbar, best, temp, rand, evals
REAL(DP), DIMENSION(1000000) :: geni,genavg,genmax
! REAL(DP) :: TIME0,cpu0,cpu1,cpu
! REAL(DP), DIMENSION(2) :: tarray
INTEGER(K4), DIMENSION(nchrmax) :: ibest
!$ INTEGER(K4) :: omp_get_max_threads, omp_get_thread_num
!
numthreads = 1
!$     numthreads=omp_get_max_threads()
!$OMP PARALLEL DO ORDERED
! create an array of locks
!$     DO i = 1,indmax
!$        CALL omp_init_lock(lck(i))
!$     ENDDO
!$OMP END PARALLEL DO 
!
! CALL etime(tarray)
! WRITE(6,*) tarray(1),tarray(2)
! cpu0=tarray(1)
!
! Call the input subroutine.
! TIME0=SECNDS(0.0)
CALL input
!
! Perform necessary initialization and read the ga.restart file.
CALL initial(istart,npossum,ig2sum)
!
!!!!!!! Main generational processing loop. !!!!!!!!!!!
kount=0
DO i=istart,maxgen+istart-1
!! WRITE (6,1111) i
!! WRITE (24,1111) i
!! WRITE(24,1050)
!
!  Evaluate the population, assign fitness, establish the best
!  individual, and write output information.
   CALL evalout(iskip,iend,ibest,fbar,best)
   geni(i)=FLOAT(i)
   genavg(i)=fbar
   genmax(i)=best
   IF(npopsiz.EQ.1 .OR. iskip.NE.0) THEN
     CLOSE(24)
     STOP
   END IF
!
!  Implement "niching".
   IF (iniche.NE.0) CALL niche
!
!  Enter selection, crossover and mutation loop.
   ncross=0
   ipick = 1
   offset = npopsiz/4  ! npopsiz should be a multiple of 4
! following loop is correct for nchild=1.
   DO k = 1,4
     CALL shuffle
! generate one-forth of the children population now
     childcnt = childcnt + offset
     crosovrcnt = crosovrcnt + offset*nchrome
     selectcnt = selectcnt + offset*nchrome
     !$OMP PARALLEL PRIVATE (j,rand,ifirst,isecond,ithird,iforth,mate1,mate2,my_cpu_id)
       my_cpu_id = 1
!$     my_cpu_id = omp_get_thread_num() + 1
     !$OMP DO SCHEDULE(GUIDED)
     DO j=ipick,ipick+offset-1,nchild
	! Perform selection.
	! CALL selectn(k,j,mate1,mate2)
	! If tournament selection is chosen (i.e. itourny=1), then
	! implement "tournament" selection for selection of new population.
	IF (itourny.EQ.1) THEN
	   ifirst=(4*(j-1)+1)-(npopsiz*(k-1))
	   isecond=ifirst+1
	   ithird=ifirst+2
	   iforth=ifirst+3
	   IF (fitness(ifirst).GT.fitness(isecond)) THEN
	      mate1=ifirst
	   ELSE
	      mate1=isecond
	   END IF
	   IF (fitness(ithird).GT.fitness(iforth)) THEN
	      mate2=ithird
	   ELSE
	      mate2=iforth
	   END IF
	   DO n=1,nchrome
	      ichild(n,j)=iparent(n,mate1)
	      IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate2)
	   END DO
	END IF
	! Now perform crossover between the randomly selected pair.
	! CALL crosovr(ncross,j,mate1,mate2,my_cpu_id)
	IF (iunifrm.EQ.0) THEN
	!  Single-point crossover at a random chromosome point.
	   CALL ran3(1,rand,my_cpu_id,0)
	   IF (rand.GT.pcross) GOTO 69
	   ncross=ncross+1
	   CALL ran3(1,rand,my_cpu_id,0)
	   icross=2+DINT(DBLE(nchrome-1)*rand)
	   DO n=icross,nchrome
	      ichild(n,j)=iparent(n,mate2)
	      IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate1)
	   END DO
	ELSE
	!  Perform uniform crossover between the randomly selected pair.
	   DO n=1,nchrome
	      CALL ran3(1,rand,my_cpu_id,0)
	      IF (rand.LE.pcross) THEN
		 ichild(n,j)=iparent(n,mate2)
		 IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate1)
	      END IF
	   END DO
	END IF
 69     CONTINUE
     END DO
     !$OMP END DO
     !$OMP END PARALLEL
     ipick = ipick + offset
   END DO
!
  maincnt = maincnt + 1
! Now perform random mutations.  If running micro-GA, skip mutation.
!
  IF (microga.EQ.0) CALL mutate
!
! Write child array back into parent array for new generation.  Check
! to see IF the best parent was replicated.
!
  CALL newgen(ielite,npossum,ig2sum,ibest)
!
!  Implement micro-GA IF enabled.
!
  IF (microga.NE.0) CALL gamicro(i,npossum,ig2sum,ibest)
!
!  Write to restart file.
!         CALL restart(i,istart,kount)
END DO
!!!!!!! End of main generational processing loop. !!!!!!
! 999 CONTINUE 
WRITE(24,3000)
DO i=1,maxgen
   evals=FLOAT(npopsiz)*geni(i)
   WRITE(24,3100) geni(i),evals,genavg(i),genmax(i)
END DO
! write validation counters
WRITE(26,5100) 
WRITE(26,5200) maincnt
WRITE(26,5200) evaloutcnt
WRITE(26,5200) childcnt
WRITE(26,5200) crosovrcnt
WRITE(26,5200) selectcnt
WRITE(26,5200) shufflecnt
WRITE(26,5200) newgencnt
WRITE(26,5200) mutatecnt
WRITE(26,5300) genmax(maxgen)
! CALL etime(tarray)
! WRITE(6,*) tarray(1),tarray(2)
! cpu1=tarray(1)
! cpu=(cpu1-cpu0)
! WRITE(6,1400) cpu,cpu/60.0
! WRITE(24,1400) cpu,cpu/60.0
CLOSE (24)
!! 1050 FORMAT(1X,' #      Binary Code',16x,'Param1  Param2  Fitness')
!! 1111 FORMAT(//'#################  Generation',I5,'  #################')
!! 1225 FORMAT(/'  Number of Crossovers      =',I5)
! 1400 FORMAT(2X,'CPU time for all generations=',E12.6,' sec'/  &
!           2X,'                             ',E12.6,' min')
3000 FORMAT(2X,'Generation   Evaluations   Avg.Fitness   Best Fitness')
3100 FORMAT(2X,4(E10.4,3X))
!
5100 FORMAT(2X,'Gafort Validation Counters and Best Fitness')
5200 FORMAT(2X,I20)
5300 FORMAT(2X,F10.9)
!
! destroy an array of locks
!$ DO i = 1,indmax
!$    CALL omp_destroy_lock(lck(i))
!$ ENDDO 
IF (lenran > 0) THEN
  DEALLOCATE(procstates)
  lenran = 0
END IF
CLOSE(26)
!
END PROGRAM gafortran
!***********************************************************************
SUBROUTINE input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine inputs information from the ga.inp (gafort.in) file. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,nposibl,nichflg,&
                           parmax,parmin,pardel,npopsiz,nowrite,      &
                           nparam,nchrome,pcross,pmutate,pcreep,      &
                           maxgen,idum,irestrt,itourny,ielite,icreep, &
                           iunifrm,iniche,iskip,iend,nchild,microga,  &
                           kountmx
IMPLICIT NONE
INTEGER(K4) :: i
NAMELIST /ga/ irestrt,npopsiz,pmutate,maxgen,idum,pcross,        &
              itourny,ielite,icreep,pcreep,iunifrm,iniche,       &
              iskip,iend,nchild,nparam,parmin,parmax,nposibl,    &
              nowrite,nichflg,microga,kountmx
!
kountmx=5
irestrt=0
itourny=0
ielite=0
iunifrm=0
iniche=0
iskip=0
iend=0
nchild=1
nichflg(1:nparam)=1
microga=0
!
OPEN (UNIT=24, FILE='gafort.out', STATUS='UNKNOWN')
REWIND 24
OPEN (UNIT=23, FILE='gafort.in', STATUS='UNKNOWN')
READ (23, NML = ga)
CLOSE (23)
OPEN (UNIT=26, FILE='gafort.vld', STATUS='UNKNOWN')
REWIND 26
itourny=1
irestrt=0
kountmx=maxgen
! IF (itourny.EQ.0) nchild=2
! Check for array sizing errors.
IF (npopsiz.GT.indmax) THEN
   WRITE(6,1600) npopsiz
   WRITE(24,1600) npopsiz
   CLOSE(24)
   STOP
END IF
IF (nparam.GT.nparmax) THEN
   WRITE(6,1700) nparam
   WRITE(24,1700) nparam
   CLOSE(24)
   STOP
END IF
! If using the microga option, reset some input variables
IF (microga.NE.0) THEN
   pmutate=0.0d0
   pcreep=0.0d0
   itourny=1
   ielite=1
   iniche=0
   nchild=1
   IF (iunifrm.EQ.0) THEN
      pcross=1.0d0
   ELSE
      pcross=0.5d0
   END IF
END IF
!
1600 FORMAT(1X,'ERROR: npopsiz > indmax.  Set indmax = ',I6)
1700 FORMAT(1X,'ERROR: nparam > nparmax.  Set nparmax = ',I6)
END SUBROUTINE input 
!***********************************************************************
SUBROUTINE initial(istart,npossum,ig2sum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets up the program by generating the g0, g1 and    !
! ig2 arrays, and counting the number of chromosomes required for the !
! specified input.  The subroutine also initializes the random number !
! generator, parent and iparent arrays (reads the ga.restart file).   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,parent,iparent,nposibl, &
                          g0,g1,ig2,parmax,parmin,pardel,npopsiz,nowrite,&
                          nparam,nchrome,pcross,pmutate,pcreep,maxgen,   &
                          idum,irestrt,itourny,ielite,icreep,iunifrm,    &
                          iniche,iskip,iend,nchild,microga,kountmx,     &
                          ichild,child,fitness,numthreads
IMPLICIT NONE
INTEGER(K4), INTENT(OUT) :: istart, npossum, ig2sum
INTEGER(K4) :: i, j, n2j
INTEGER(K4), DIMENSION(numthreads) :: iduma
REAL(DP) :: rand
INTEGER(K4), DIMENSION(nparam) :: ig2sum1
!
g0(1:nparam)=parmin(1:nparam)
pardel(1:nparam)=parmax(1:nparam)-parmin(1:nparam)
g1(1:nparam)=pardel(1:nparam)/DBLE(nposibl(1:nparam)-1)
DO i=1,nparam
   DO j=1,30
      n2j=2**j
      IF (n2j.GE.nposibl(i)) THEN
         ig2(i)=j
         GOTO 8
      END IF
      IF (j.GE.30) THEN
         WRITE(6,2000)
         WRITE(24,2000)
         CLOSE(24)
         STOP
      END IF
   END DO
8  CONTINUE
END DO
! Count the total number of chromosomes (bits) required
nchrome=SUM(ig2(1:nparam))
npossum=SUM(nposibl(1:nparam))
ig2sum1(1:nparam)=2**ig2(1:nparam)
ig2sum=SUM(ig2sum1(1:nparam))
IF (nchrome.GT.nchrmax) THEN
   WRITE(6,1800) nchrome
   WRITE(24,1800) nchrome
   CLOSE(24)
   STOP
END IF
IF (npossum.LT.ig2sum .AND. microga.NE.0) THEN
   WRITE(6,2100)
   WRITE(24,2100)
END IF
!
!distribute arrays across memories
!$OMP PARALLEL
!$OMP DO
   DO i=1,npopsiz
      fitness(i)=0
      DO j=1,nchrome
         iparent(j,i)=0
         ichild(j,i)=0
      ENDDO
      DO j=1,nparam
         parent(j,i)=0
         child(j,i)=0
      END DO
   END DO    
!$OMP END DO
!$OMP END PARALLEL
!
! Initialize random number generator
CALL ran3(idum,rand,1,numthreads)
!
IF (irestrt.EQ.0) THEN
!  Initialize the random distribution of parameters in the individual
!  parents when irestrt=0.
   istart=1
   DO i=1,npopsiz
      DO j=1,nchrome
         CALL ran3(1,rand,1,0)
         iparent(j,i)=1
         IF(rand.LT.0.5d0) iparent(j,i)=0
      END DO
   END DO
   IF (npossum.LT.ig2sum) CALL possibl(parent,iparent)
END IF
IF(irestrt.NE.0) CALL ran3(idum-istart,rand,1,numthreads)
1800 FORMAT(1X,'ERROR: nchrome > nchrmax.  Set nchrmax = ',I6)
2000 FORMAT(1X,'ERROR: You have a parameter with a number of '/           &
            1X,'   possibilities > 2**30!  If you really desire this,'/   &
            1X,'   change the DO loop 7 statement and recompile.'//       &
            1X,'   You may also need to alter the code to work with'/     &
            1X,'   REAL numbers rather than INTEGER numbers; Fortran'/    &
            1X,'   does not like to compute 2**j when j>30.')             
2100 FORMAT(1X,'WARNING: for some cases, a considerable performance'/     &
            1X,'   reduction has been observed when running a non-'/      &
            1X,'   optimal number of bits with the micro-GA.'/            &
            1X,'   If possible, use values for nposibl of 2**n,'/         &
            1X,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
END SUBROUTINE initial
!***********************************************************************
SUBROUTINE evalout(iskip,iend,ibest,fbar,best)
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,parent,iparent,   &
                          g0,g1,ig2,fitness,npopsiz,nowrite,nparam, &
                          nchrome,numthreads,evaloutcnt
IMPLICIT NONE
SAVE
INTEGER(K4), INTENT(IN) :: iskip,iend
INTEGER(K4), INTENT(OUT), DIMENSION(nchrmax) :: ibest
REAL(DP), INTENT(OUT) :: fbar,best
REAL(DP), DIMENSION(nparmax) :: paramav,paramsm
INTEGER(K4) :: j, n, k, jstart, jend, jbest, my_cpu_id
REAL(DP) :: fitsum, pi
REAL(DP), DIMENSION(nparmax,numthreads) :: paramsm1
REAL(DP), DIMENSION(nparam) :: f1,f2,funcval1
INTEGER(k4) :: l, iparam, m, i, nvalley
!$ INTEGER(K4) :: omp_get_thread_num
!
fitsum=0.0d0
best=-1.0d10
paramsm(1:nparam)=0.0D0
jstart=1
jend=npopsiz
IF(iskip.NE.0) jstart=iskip
IF(iend.NE.0) jend=iend
paramsm1(1:nparmax,1:numthreads) = 0.0
!$OMP PARALLEL PRIVATE(j,my_cpu_id,n,l,k,iparam,m,nvalley,pi,f1,f2,funcval1) REDUCTION(+:fitsum)
my_cpu_id = 1
!$ my_cpu_id = omp_get_thread_num() + 1
!$OMP DO
DO j=jstart,jend
!   CALL decode(j,parent,iparent)
    l=1
    DO k=1,nparam
       iparam=0
       m=l
       DO i=m,m+ig2(k)-1
          l=l+1
          iparam=iparam+iparent(i,j)*(2**(m+ig2(k)-1-i))
       END DO
       parent(k,j)=g0(k)+g1(k)*DBLE(iparam)
    END DO
!! IF(iskip.NE.0 .AND. iend.NE.0 .AND. iskip.EQ.iend) &
!     WRITE(6,1075) j,(iparent(k,j),k=1,nchrome),    &
!!                 (parent(kk,j),kk=1,nparam),0.0
!
!  Call function evaluator, write out individual and fitness, and add
!  to the summation for later averaging.
!  CALL func(j,funcval)
   nvalley=4
   pi=4.0d0*datan(1.d0)
   f1(1:nparam)=(SIN(5.1d0*pi*parent(1:nparam,j) + 0.5d0))**nvalley
   f2(1:nparam)=EXP(-4.0d0*log(2.0d0)*((parent(1:nparam,j)-0.0667d0)**2)/0.64d0)
   funcval1(1:nparam)=f1(1:nparam)*f2(1:nparam)
   fitness(j)=PRODUCT(funcval1(1:nparam))
!!         write(24,1075) j,(iparent(k,j),k=1,nchrome),      &
!!                  (parent(kk,j),kk=1,nparam),fitness(j)
   fitsum=fitsum+fitness(j)
   paramsm1(1:nparam,my_cpu_id)=paramsm1(1:nparam,my_cpu_id)+parent(1:nparam,j)
END DO
!$OMP END DO
!$OMP END PARALLEL
!
! Check to see if fitness of individual j is the best fitness.
DO j=jstart,jend
   IF (fitness(j).GT.best) THEN 
      best=fitness(j)
      jbest=j
      ibest(1:nchrome)=iparent(1:nchrome,j)
   END IF
END DO
DO j=1,numthreads
   paramsm(1:nparam)=paramsm(1:nparam)+paramsm1(1:nparam,j)
END DO
evaloutcnt = evaloutcnt + (jend-jstart)
fbar=fitsum/DBLE(npopsiz)
paramav(1:nparam)=paramsm(1:nparam)/DBLE(npopsiz)
!
!  Write output information
!! IF (npopsiz.EQ.1) THEN 
!!     WRITE(24,1075) 1,(iparent(k,1),k=1,nchrome),  &
!!           (parent(k,1),k=1,nparam),fitness(1)
!!     WRITE(24,*) ' Average Values:'
!!     WRITE(24,1275) (parent(k,1),k=1,nparam),fbar
!! ELSE 
!!     WRITE(24,1275) (paramav(k),k=1,nparam),fbar
!! END IF
!! WRITE(6,1100) fbar
!! WRITE(24,1100) fbar
!! WRITE(6,1200) best
!! WRITE(24,1200) best
!
!! 1075 FORMAT(I3,1X,30X1,2(2X,F6.3),2X,F6.3)
!! 1100 FORMAT(1X,'Average Function Value of Generation=',F6.3)
!! 1200 FORMAT(1X,'Maximum Function Value              =',F6.3/)
!! 1275 FORMAT(/' Average Values:',18X,2(2X,F6.3),2X,F6.3/)
END SUBROUTINE evalout
!*************************************************************************
SUBROUTINE niche
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Implement "niching" through Goldberg's multidimensional phenotypic    !
!  sharing scheme with a triangular sharing function.  To find the       !
!  multidimensional distance from the best individual, normalize all     !
!  parameter differences.                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Variable definitions:                                                !
!                                                                        !
!  alpha   = power law exponent for sharing function; typically = 1.0    !
!  del     = normalized multidimensional distance between ii and all     !
!            other members of the population                             !
!            (equals the square root of del2)                            !
!  del2    = sum of the squares of the normalized multidimensional       !
!            distance between member ii and all other members of         !
!            the population                                              !
!  nniche  = number of niched parameters                                 !
!  sigshar = normalized distance to be compared with del; in some sense, !
!            1/sigshar can be viewed as the number of regions over which !
!            the sharing function should focus, e.g. with sigshar=0.1,   !
!            the sharing function will try to clump in ten distinct      !
!            regions of the phase space.  A value of sigshar on the      !
!            order of 0.1 seems to work best.                            !
!  share   = sharing function between individual ii and j                !
!  sumshar = sum of the sharing functions for individual ii              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,parent,iparent,fitness, &
                           nposibl,nichflg,parmax,parmin,pardel,npopsiz,  &
                           nowrite,nparam,nchrome
IMPLICIT NONE
SAVE
INTEGER(K4) :: nniche, ii, j, k, jj
! REAL(DP) :: alpha
REAL(DP) :: sigshar, sumshar, del2, del, share
! alpha=1.0
sigshar=0.1d0
nniche=SUM(nichflg(1:nparam))
IF (nniche.EQ.0) THEN
   WRITE(6,1900)
   WRITE(24,1900)
   CLOSE(24)
   STOP
END IF
!$OMP PARALLEL PRIVATE(ii,sumshar,j,del2,k,del,share)
!$OMP DO
DO ii=1,npopsiz
   sumshar=0.0d0
   DO j=1,npopsiz
      del2=0.0d0
      DO k=1,nparam
         IF (nichflg(k).NE.0) THEN
            del2=del2+((parent(k,j)-parent(k,ii))/pardel(k))**2
         END IF
      END DO
      del=(DSQRT(del2))/DBLE(nniche)
      IF (del.LT.sigshar) THEN
!        share=1.0-((del/sigshar)**alpha)
         share=1.0d0-(del/sigshar)
      ELSE
         share=0.0d0
      END IF
      sumshar=sumshar+share/DBLE(npopsiz)
   END DO
   IF (sumshar.NE.0.0d0) fitness(ii)=fitness(ii)/sumshar
END DO
!$OMP END DO
!$OMP END PARALLEL
!
1900 FORMAT(1X,'ERROR: iniche=1 and all values in nichflg array = 0'/   &
            1X,'       Do you want to niche or not?')
END SUBROUTINE niche
!*************************************************************************
SUBROUTINE selectn(k,j,mate1,mate2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for selection operator.  Presently, tournament selection !
! is the only option available.                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,fitness,iparent,ichild,npopsiz,nchrome, &
                           itourny,nchild
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: j,k 
INTEGER(K4), INTENT(OUT) :: mate1,mate2 
INTEGER(K4) :: n,ifirst,isecond,ithird,iforth
! If tournament selection is chosen (i.e. itourny=1), then
! implement "tournament" selection for selection of new population.
IF (itourny.EQ.1) THEN
   ifirst=(4*(j-1)+1)-(npopsiz*(k-1))
   isecond=ifirst+1
   ithird=ifirst+2
   iforth=ifirst+3
! select first mate
   IF (fitness(ifirst).GT.fitness(isecond)) THEN
      mate1=ifirst
   ELSE
      mate1=isecond
   END IF
! select second mate
   IF (fitness(ithird).GT.fitness(iforth)) THEN
      mate2=ithird
   ELSE
      mate2=iforth
   END IF
!  WRITE(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
   DO n=1,nchrome
      ichild(n,j)=iparent(n,mate1)
      IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate2)
   END DO
END IF
END SUBROUTINE selectn
!*************************************************************************
SUBROUTINE crosovr(ncross,j,mate1,mate2,my_cpu_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for crossover between the randomly selected pair. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,iparent,ichild,nchrome,pcross,iunifrm,nchild
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: j, mate1, mate2,my_cpu_id
INTEGER(K4), INTENT(OUT) :: ncross
INTEGER(K4) :: icross, n
REAL(DP) :: rand
!
IF (iunifrm.EQ.0) THEN
!  Single-point crossover at a random chromosome point.
   CALL ran3(1,rand,my_cpu_id,0)
   IF (rand.GT.pcross) GOTO 69
   ncross=ncross+1
   CALL ran3(1,rand,my_cpu_id,0)
   icross=2+DINT(DBLE(nchrome-1)*rand)
   DO n=icross,nchrome
      ichild(n,j)=iparent(n,mate2)
      IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate1)
   END DO
ELSE
!  Perform uniform crossover between the randomly selected pair.
   DO n=1,nchrome
      CALL ran3(1,rand,my_cpu_id,0)
      IF (rand.LE.pcross) THEN
         ncross=ncross+1
         ichild(n,j)=iparent(n,mate2)
         IF(nchild.EQ.2) ichild(j+1,n)=iparent(n,mate1)
      END IF
   END DO
END IF
69   CONTINUE
END SUBROUTINE crosovr
!*************************************************************************
SUBROUTINE mutate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs mutations on the children generation.        !
! Perform random jump mutation if a random number is less than pmutate. !
! Perform random creep mutation if a different random number is less    !
! than pcreep.                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           npopsiz,nowrite,nparam,nchrome,g0,g1,ig2, &
                           parmax,parmin,pardel,nposibl,child,ichild, &
                           pcross,pmutate,pcreep,maxgen,idum,irestrt, &
                           itourny,ielite,icreep,iunifrm,iniche,   &
                           iskip,iend,nchild,microga,kountmx,mutatecnt
IMPLICIT NONE
INTEGER(K4) :: nmutate, ncreep
INTEGER(K4) :: j, k, my_cpu_id
!$ INTEGER(K4) :: omp_get_thread_num
INTEGER(K4) :: l,m,i,kk,iparam
REAL(DP) :: rand, creep, temp(nchrmax)
nmutate=0
ncreep=0
!$OMP PARALLEL PRIVATE(j,k,rand,my_cpu_id,temp,nmutate)
      my_cpu_id = 1
!$    my_cpu_id = omp_get_thread_num() + 1
!$OMP DO  SCHEDULE(GUIDED)
DO j=1,npopsiz
   DO k=1,nchrome
!     Jump mutation
      CALL ran3(1,rand,my_cpu_id,0)
      IF (rand.LE.pmutate) THEN
         nmutate=nmutate+1
         IF (ichild(k,j).EQ.0) THEN
            ichild(k,j)=1
         ELSE
            ichild(k,j)=0
         END IF
         IF (nowrite.EQ.0) WRITE(6,1300) j,k
         IF (nowrite.EQ.0) WRITE(24,1300) j,k
      END IF
   END DO
!  Creep mutation (one discrete position away).
   IF (icreep.NE.0) THEN
      DO k=1,nparam
         CALL ran3(1,rand,my_cpu_id,0)
         IF (rand.LE.pcreep) THEN
            CALL decode(j,child,ichild)
            ncreep=ncreep+1
            creep=1.0d0
            CALL ran3(1,rand,my_cpu_id,0)
            IF (rand.LT.0.5d0) creep=-1.0d0
            child(k,j)=child(k,j)+g1(k)*creep
            IF (child(k,j).GT.parmax(k)) THEN
               child(k,j)=parmax(k)-1.0d0*g1(k)
            ELSE IF (child(k,j).LT.parmin(k)) THEN
               child(k,j)=parmin(k)+1.0d0*g1(k)
            END IF
            CALL code(j,k,child,ichild)
            IF (nowrite.EQ.0) WRITE(6,1350) j,k
            IF (nowrite.EQ.0) WRITE(24,1350) j,k
         END IF
      END DO  
   END IF
END DO
!$OMP END DO
!$OMP END PARALLEL
mutatecnt = mutatecnt + npopsiz
!WRITE(6,1250) nmutate,ncreep
!WRITE(26,1250) nmutate,ncreep
!
!1250 FORMAT(/'  Number of Jump Mutations  =',I5/  &
!           '  Number of Creep Mutations =',I5)
1300 FORMAT('*** Jump mutation performed on individual  ',I4,  &
            ', chromosome ',I3,' ***')
1350 FORMAT('*** Creep mutation performed on individual ',I4,  &
            ', parameter  ',I3,' ***')
!
END SUBROUTINE mutate
!*************************************************************************
SUBROUTINE newgen(ielite,npossum,ig2sum,ibest)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write child array back into parent array for new generation.  Check !
! to see if the best parent was replicated; if not, and if ielite=1,  !
! then reproduce the best parent into a random slot.                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           npopsiz,nowrite,nparam,nchrome, &
                           parent,iparent,child,ichild,numthreads, &
			   newgencnt
IMPLICIT NONE
SAVE
INTEGER(K4), INTENT(IN) :: ielite, npossum, ig2sum
INTEGER(K4), DIMENSION(nchrmax) :: ibest
INTEGER(K4) :: j, n, kelite, jelite, irand, my_cpu_id
!$ INTEGER(K4) :: omp_get_thread_num
INTEGER(K4), DIMENSION(numthreads) :: kelite1
REAL(DP) :: rand
!
IF (npossum.LT.ig2sum) CALL possibl(child,ichild)
!$OMP PARALLEL PRIVATE(j,n,my_cpu_id)
my_cpu_id = 1
!$ my_cpu_id = omp_get_thread_num() + 1
kelite1(my_cpu_id) = 0
!$OMP DO LASTPRIVATE(jelite,kelite)
  DO j=1,npopsiz
     jelite=0
     DO n=1,nchrome
        iparent(n,j)=ichild(n,j)
        IF (iparent(n,j).EQ.ibest(n)) jelite=jelite+1
        IF (jelite.EQ.nchrome) kelite1(my_cpu_id)=1
     END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
kelite=SUM(kelite1) ! If any elements are 1, kelite>0, else kelite=0.
newgencnt = newgencnt + npopsiz*nchrome
IF (ielite.NE.0 .AND. kelite.EQ.0) THEN
   CALL ran3(1,rand,1,0)
   irand=1d0+DINT(DBLE(npopsiz)*rand)
   iparent(1:nchrome,irand)=ibest(1:nchrome)
!  WRITE(24,1260) irand
END IF
!
! 1260 FORMAT('  Elitist Reproduction on Individual ',I4)
!
END SUBROUTINE newgen
!*************************************************************************
SUBROUTINE gamicro(i,npossum,ig2sum,ibest)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Micro-GA implementation subroutine !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           npopsiz,nowrite,nparam,nchrome, &
                           parent,iparent
IMPLICIT NONE
SAVE
INTEGER(K4), INTENT(IN) :: i, npossum, ig2sum
INTEGER(K4), DIMENSION(nchrmax) :: ibest
INTEGER(K4) :: icount, j, n
REAL(DP) :: diffrac, rand
! First, check for convergence of micro population.
! If converged, start a new generation with best individual and fill
! the remainder of the population with new randomly generated parents.
!
! Count number of different bits from best member in micro-population
icount=0
DO j=1,npopsiz
   DO n=1,nchrome
      IF(iparent(n,j).NE.ibest(n)) icount=icount+1
   END DO
END DO
!  If icount less than 5% of number of bits, then consider population
!  to be converged.  Restart with best individual and random others.
diffrac=DBLE(icount)/DBLE((npopsiz-1)*nchrome)
IF (diffrac.LT.0.05d0) THEN
   iparent(1:nchrome,1)=ibest(1:nchrome)
   DO j=2,npopsiz
      DO n=1,nchrome
         CALL ran3(1,rand,1,0)
         iparent(n,j)=1
         IF(rand.LT.0.5d0) iparent(n,j)=0
      END DO 
   END DO 
   IF (npossum.LT.ig2sum) CALL possibl(parent,iparent)
!  WRITE(6,1375) i
!  WRITE(24,1375) i
END IF
!
! 1375 FORMAT(//'%%%%%%%  Restart micro-population at generation',    &
!            I5,'  %%%%%%%')
END SUBROUTINE gamicro 
!*************************************************************************
SUBROUTINE shuffle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine shuffles the parent array and its corresponding fitness !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           npopsiz,nowrite,nparam,nchrome, &
                           parent,iparent,fitness,lck,shufflecnt
IMPLICIT NONE
INTEGER(K4) :: j, itemp(nchrome), iother
!$ INTEGER(K4) :: omp_get_thread_num
INTEGER(K4) :: my_cpu_id
REAL(DP) :: rand, temp
!$OMP PARALLEL PRIVATE(rand, iother, itemp, temp, my_cpu_id) 
my_cpu_id = 1
!$ my_cpu_id = omp_get_thread_num() + 1
!$OMP DO
DO j=1,npopsiz-1
   CALL ran3(1,rand,my_cpu_id,0)
   iother=j+1+DINT(DBLE(npopsiz-j)*rand)
!$   IF (j < iother) THEN
!$      CALL omp_set_lock(lck(j))
!$      CALL omp_set_lock(lck(iother))
!$   ELSE
!$      CALL omp_set_lock(lck(iother))
!$      CALL omp_set_lock(lck(j))
!$   END IF
   itemp(1:nchrome)=iparent(1:nchrome,iother)
   iparent(1:nchrome,iother)=iparent(1:nchrome,j)
   iparent(1:nchrome,j)=itemp(1:nchrome)
   temp=fitness(iother)
   fitness(iother)=fitness(j)
   fitness(j)=temp
!$   IF (j < iother) THEN
!$      CALL omp_unset_lock(lck(iother))
!$      CALL omp_unset_lock(lck(j))
!$   ELSE
!$      CALL omp_unset_lock(lck(j))
!$      CALL omp_unset_lock(lck(iother))
!$   END IF
END DO
!$OMP END DO
!$OMP END PARALLEL
shufflecnt = shufflecnt + npopsiz - 1
END SUBROUTINE shuffle 
!*************************************************************************
SUBROUTINE decode(i,array,iarray)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine decodes a binary string to a real number.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           nparam,nchrome,g0,g1,ig2
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: i
INTEGER(K4), DIMENSION(nchrmax,indmax) :: iarray
REAL(DP), INTENT(OUT), DIMENSION(nparmax,indmax) :: array
INTEGER(k4) :: l, k, iparam, j, m
l=1
DO k=1,nparam
   iparam=0
   m=l
   DO j=m,m+ig2(k)-1
      l=l+1
      iparam=iparam+iarray(j,i)*(2**(m+ig2(k)-1-j))
   END DO
   array(k,i)=g0(k)+g1(k)*DBLE(iparam)
END DO
END SUBROUTINE decode
!*************************************************************************
SUBROUTINE code(j,k,array,iarray)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine codes a parameter into a binary string.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax, &
                           nparam,nchrome,g0,g1,ig2
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: j, k
INTEGER(K4), INTENT(OUT), DIMENSION(nchrmax,indmax) :: iarray
REAL(DP), DIMENSION(nparmax,indmax) :: array
INTEGER(K4) :: i, iparam, istart, m
! First, establish the beginning location of the parameter string of
! interest.
istart=1+SUM(ig2(1:k-1))
!  Find the equivalent coded parameter value, and back out the binary
!  string by factors of two.
m=ig2(k)-1
IF (g1(k).EQ.0.0d0) RETURN
iparam=nint((array(k,j)-g0(k))/g1(k))
DO i=istart,istart+ig2(k)-1
   iarray(i,j)=0
   IF ((iparam+1).GT.(2**m)) THEN
      iarray(i,j)=1
      iparam=iparam-2**m
   END IF
   m=m-1
END DO
!  WRITE(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
END SUBROUTINE code
!*************************************************************************
SUBROUTINE possibl(array,iarray)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This subroutine determines whether or not all parameters are within !
!  the specified range of possibility.  If not, the parameter is       !
!  randomly reassigned within the range.  This subroutine is only      !
!  necessary when the number of possibilities per parameter is not     !
!  optimized to be 2**n, i.e. if npossum < ig2sum.                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,npopsiz,nowrite,  &
                           nparam,nchrome,g0,g1,ig2,  &
                           parmax,parmin,pardel,nposibl
IMPLICIT NONE
REAL(DP), INTENT(OUT), DIMENSION(nparmax,indmax) :: array
INTEGER(K4), INTENT(INOUT), DIMENSION(nchrmax,indmax) :: iarray
INTEGER(K4) :: i, j, n2ig2j, irand, my_cpu_id
!$ INTEGER(K4) :: omp_get_thread_num
REAL(DP) :: rand
!$OMP PARALLEL PRIVATE(i,j,my_cpu_id,n2ig2j,irand,rand)
my_cpu_id = 1
!$ my_cpu_id = omp_get_thread_num() + 1
!$OMP DO SCHEDULE(GUIDED)
DO i=1,npopsiz
   CALL decode(i,array,iarray)
   DO j=1,nparam
      n2ig2j=2**ig2(j)
      IF (nposibl(j).NE.n2ig2j .AND. array(j,i).GT.parmax(j)) THEN
         CALL ran3(1,rand,1,0)
         irand=DINT(DBLE(nposibl(j))*rand)
         array(j,i)=g0(j)+DBLE(irand)*g1(j)
         CALL code(i,j,array,iarray)
         IF (nowrite.EQ.0) WRITE(6,1000) i,j
         IF (nowrite.EQ.0) WRITE(24,1000) i,j
      END IF
   END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
1000 FORMAT('*** Parameter adjustment to individual     ',I4,   &
            ', parameter  ',I3,' ***')
END SUBROUTINE possibl
!*************************************************************************
! SUBROUTINE restart(i,istart,kount)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine writes restart information to the ga.restart file.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,npopsiz,nowrite, &
!                             nparam,nchrome,parent,iparent,pcross,  &
!                             pmutate,pcreep,maxgen,idum,irestrt,itourny, &
!                             ielite,icreep,iunifrm,iniche,iskip,iend, &
!                             nchild,microga,kountmx
! IMPLICIT NONE
! SAVE
! INTEGER(K4), INTENT(IN) :: i, istart
! INTEGER(K4), INTENT(OUT) :: kount
! kount=kount+1
! IF (i.EQ.maxgen+istart-1 .OR. kount.EQ.kountmx) THEN
!    OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
!    REWIND 25
!    WRITE(25,*) i+1,npopsiz
!    DO j=1,npopsiz
!       WRITE(25,1500) j,(iparent(l,j),l=1,nchrome)
!    END DO
!    CLOSE (25)
!    kount=0
! END IF
!
! 1500 FORMAT(I5,3X,30I2)
! END SUBROUTINE restart
!*************************************************************************
SUBROUTINE ran_init(idum,procs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the state of the random number generator per processor.   !
! To avoid false sharing and increase spatial locality, we allocate    !
! one structure ranstate per processor and maintain pointers to it.    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP
USE ran_state
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: procs,idum
INTEGER(K4) :: i,ii,k,j
REAL(DP) :: mk
!
! See if we need to deallocate space.
!
IF (lenran > 0) THEN
   DEALLOCATE(procstates)
   lenran=0
END IF
!
! Allocate memory for the state variables.
!
lenran=procs
ALLOCATE(procstates(lenran))    ! Allocate an array of pointers to ranstate per proc.
DO j=1,lenran
   procstates(j)%iff=0
   procstates(j)%iduma=idum*j  ! Generate seed value for each processor.
   IF (procstates(j)%iff.EQ.0) THEN
      procstates(j)%iff=1
      procstates(j)%mj=mseed-DBLE(IABS(procstates(j)%iduma))
      procstates(j)%mj=DMOD(procstates(j)%mj,mbig)
      procstates(j)%ma(madim)=procstates(j)%mj
      mk=1
      DO i=1,madim-1
         ii=MOD(21*i,madim)
         procstates(j)%ma(ii)=mk
         mk=procstates(j)%mj-mk
         IF(mk.LT.mz) mk=mk+mbig
         procstates(j)%mj=procstates(j)%ma(ii)
      END DO
      DO k=1,4
         DO i=1,madim
            procstates(j)%ma(i)=procstates(j)%ma(i)-procstates(j)%ma(1+MOD(i+30,madim))
            IF(procstates(j)%ma(i).LT.mz) procstates(j)%ma(i)=procstates(j)%ma(i)+mbig
         END DO
      END DO
      procstates(j)%inext=0
      procstates(j)%inextp=31
      procstates(j)%iduma=1
   END IF
END DO
END SUBROUTINE ran_init
!*************************************************************************
SUBROUTINE ran3(idum,rand,j,numprocs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to  !
!  any negative value to initialize or reinitialize the sequence.      !
!  This function is taken from W.H. Press', "Numerical Recipes" p. 199.!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE gafort_globals, ONLY : K4,DP
USE ran_state
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: j,numprocs,idum
REAL(DP), INTENT(OUT) :: rand
!
IF (idum.LT.0) CALL ran_init(idum,numprocs)
procstates(j)%inext=procstates(j)%inext+1
IF(procstates(j)%inext.EQ.56) procstates(j)%inext=1
procstates(j)%inextp=procstates(j)%inextp+1
IF(procstates(j)%inextp.EQ.56) procstates(j)%inextp=1
procstates(j)%mj=procstates(j)%ma(procstates(j)%inext) & 
                       - procstates(j)%ma(procstates(j)%inextp)
IF(procstates(j)%mj.LT.mz) procstates(j)%mj=procstates(j)%mj+mbig
procstates(j)%ma(procstates(j)%inext)=procstates(j)%mj
rand=procstates(j)%mj*fac
END SUBROUTINE ran3
!*************************************************************************
SUBROUTINE func(j,funcval)
USE gafort_globals, ONLY : K4,DP,indmax,nchrmax,nparmax,parent,iparent, &
                           nparam,nchrome
IMPLICIT NONE
INTEGER(K4), INTENT(IN) :: j
REAL(DP), INTENT(OUT) :: funcval
INTEGER(K4) :: nvalley, i
REAL(DP) :: pi
REAL(DP), DIMENSION(nparam) :: f1,f2,funcval1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This is an N-dimensional version of the multimodal function with  !
!  decreasing peaks used by Goldberg and Richardson (1987, see ReadMe! 
!  file for complete reference).  In N dimensions, this function has !
!  (nvalley-1)^nparam peaks, but only one global maximum.  It is a   !
!  reasonably tough problem for the GA.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nvalley=4
pi=4.0d0*datan(1.d0)
f1(1:nparam)=(SIN(5.1d0*pi*parent(1:nparam,j) + 0.5d0))**nvalley
f2(1:nparam)=EXP(-4.0d0*log(2.0d0)*((parent(1:nparam,j)-0.0667d0)**2)/0.64d0)
funcval1(1:nparam)=f1(1:nparam)*f2(1:nparam)
funcval=PRODUCT(funcval1(1:nparam))
END SUBROUTINE func
!*************************************************************************