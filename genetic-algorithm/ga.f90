!=======================================================
PROGRAM my_cga
        ! Main genetic algorithm program for Haupt and Haupt
        ! 2003 - uses High Performance Fortran
        ! Purpose: Optimizes variables based on a cost
        ! function using a genetic algorithm. Based on pseudocode in
        ! Haupt and Haupt, 1998
        !
        ! Date        Programmer     Description of Changes 
        ! =====       ============   ==========================
        ! 3 July 2003 Jaymon Knight  Code based on pseudocode
        ! 19 Nov 2003 Sue Haupt      Revised for 2nd ed of Haupt and Haupt
        ! June 2017   Nikola Mirkov  Adapt to "ordinary" fortran from HPF version
        !
        !
        !
        USE funct
        USE qsort_c_module

        IMPLICIT NONE

        ! Define parameters
        ! Define GA parameters
        ! Use these to tune the code to your problem

        INTEGER,PARAMETER::maxit=50     ! Maximum number of iterations
        INTEGER,PARAMETER::max_same=50  ! Maximum# of consecutively equal vals
        INTEGER,PARAMETER::popsize=100  ! Size of population
        INTEGER,PARAMETER::npar=2       ! Number of parameters
        REAL,PARAMETER::tol=.01         ! Percent error for stop criteria
        REAL,PARAMETER::mutrate=0.1     ! Mutation rate
        REAL,PARAMETER::selection=0.5   ! Fraction of population to keep
        REAL,PARAMETER::lo=0.!0.      ! Minimum parameter value
        REAL,PARAMETER::hi=10.!10.      ! Maximum parameter value

        ! Define variables
        INTEGER::status    ! Error flag
        INTEGER::how_big   ! Used in the RANDOM_SEED subroutine
        INTEGER::keep      ! Number kept from each generation
        INTEGER::M         ! Number of matings
        INTEGER::nmut      ! Total number of mutations
        INTEGER::iga       ! Generation counter
        INTEGER::i,j       ! Indices
        INTEGER::same      ! Counter for consecutively equal values
        INTEGER::bad_sort  ! Counts number of bad sorts from hpf grade_up - maybe not needed now?
        REAL::minc         ! Minimum cost
        REAL::temp         ! Temporary variable
        REAL::xy           ! Mix from ma and pa


        ! Define matrix variables
        INTEGER,ALLOCATABLE,DIMENSION(:)::vals       ! Contains values from the time/date call
        INTEGER,ALLOCATABLE,DIMENSION(:)::seeds      ! Vector w/ vals for RANDOM_SEED brtn
        INTEGER,ALLOCATABLE,DIMENSION(:)::ind        ! Sorted indices from cost function
        INTEGER,ALLOCATABLE,DIMENSION(:)::ind2       ! For sorting mutated population
        INTEGER,ALLOCATABLE,DIMENSION(:)::ma,pa      ! Parents (indices)


        INTEGER,ALLOCATABLE,DIMENSION(:)::xp         ! Crossover point
        INTEGER,ALLOCATABLE,DIMENSION(:)::ix         ! Index of mate #1
        INTEGER,ALLOCATABLE,DIMENSION(:)::mrow,mcol  ! Used for sorting mutations
        REAL,ALLOCATABLE,DIMENSION(:,:)::par,par2    ! Matrix of population values
        REAL,ALLOCATABLE,DIMENSION(:)::cost          ! Cost function evaluated
        REAL,ALLOCATABLE,DIMENSION(:)::odds,odds_tmp ! Involved in pairing
        REAL,ALLOCATABLE,DIMENSION(:)::pick1,pick2   ! Mates one and two
        REAL,ALLOCATABLE,DIMENSION(:)::temp_arr_1    ! Temporary 1-d array
        REAL,ALLOCATABLE,DIMENSION(:)::r             ! Mixing parameter

        ! Calculate variables
        keep=FLOOR(selection*popsize)                ! Number to keep from each generation
        M=CEILING(REAL(popsize-keep)*0.5)             ! Number of matings
        nmut=CEILING((popsize-1)*npar*mutrate)       ! Number of mutations

        !Allocate arrays (block 1)
        ALLOCATE(cost(popsize),par(popsize,npar),par2(popsize, npar),ind(popsize),odds(keep+1),odds_tmp(keep+1), &
                vals(8),ma(M),pa(M),pick1(M),pick2(M), &
                r(M),xp(M),ix(CEILING(REAL(keep)*0.5)),STAT=status)
        IF(status/=0) THEN
                WRITE(*,*) 'Error allocating arrays in main program.'
                STOP
        END IF

        !_______________________________________________________
        ! Initialize random number generator
        ! Some machines may require more care in calling the random number generator
        ! This program sets a seed randomly from computer clock

        CALL RANDOM_SEED(SIZE=how_big) ! Finds the size of array expected by subroutine
        ALLOCATE(seeds(how_big),STAT=status)
        IF(status/=0) THEN
                WRITE(*,*) '  An error occurred allocating the array ‘seeds’ in the main program.'
        END IF

        CALL DATE_AND_TIME(VALUES=vals)   !These values depend on the current time
        IF(vals(8)==0) THEN               ! We only want a non-zero value
                vals(8)=vals(5)+vals(6)+vals(7) ! Substitute in the case of zero (HH+MM+SS)
        END IF

        CALL RANDOM_SEED            ! Initializes the seed
        CALL RANDOM_SEED(GET=seeds) ! Gets the seed
        seeds=seeds*vals(8)         ! Adjusts seed so it is different each time
        CALL RANDOM_SEED(PUT=seeds) ! Seeds the random number generator

        DEALLOCATE(vals,seeds)

        !_______________________________________________________
        ! Create the initial population, evaluate costs, sort

        CALL RANDOM_NUMBER(par) ! Fills par matrix w/ random numbers 
        par=(hi-lo)*par+lo      ! between hi & lo !Normalizes values

        !_______________________________________________________
        ! Start generations

        iga=0
        minc=0.
        same=0
        bad_sort=0

        OPEN(UNIT=10,FILE='data.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=status)
        IF(status/=0) THEN
                WRITE(*,*)"Error opening file 'out.dat'."
        END IF

        DO WHILE(iga<maxit)

        iga=iga+1          ! Increment counter

        CALL ff(par,cost)  ! <<<< Calculates cost using function ff <<<<

        ! Arange (1,popsize) and store in ind
        do i=1,popsize
        ind(i) = i
        enddo

        ! Quicksort array and return sorted array and new index order relative to old one
        ! Min cost in element 1,
        ! Cost in order stored in ind
        call QsortCIndx(cost,ind)

        ! WRITE(*,*) minc,cost(1),iga

        IF(ABS((cost(1)-minc)/cost(1))<tol/100.) THEN
                same=same+1
        ELSE
                same=0
        END IF

        minc=cost(1)

        par=par(ind,:) ! Puts par in the order stored in ind.

        !_______________________________________________________
        ! Pair chromosomes and produce offspring
        odds(1)=0. !first spot is zero

        ! Fills remainder of odds matrix w/ values keep-1
        DO i=1,keep
        odds(i+1)=keep+1-i
        END DO

        !# HPF function
        !# odds(2:keep+1)=SUM_PREFIX(odds(2:keep+1))

        ! Instead of SUM_PREFIX function
        do i=2,keep+1
        odds_tmp(i) = sum(odds(2:i))
        enddo
        odds(2:keep+1)=odds_tmp(2:keep+1) ! weights chromosomes based upon position in the list

        temp=odds(keep+1)
        odds(2:keep+1)=odds(2:keep+1)/temp

        ! Probability distribution function
        CALL RANDOM_NUMBER(pick1) !mate #1
        CALL RANDOM_NUMBER(pick2) !mate #2

        ! ma and pa contain the indices of the chromosomes that will mate
        DO i=1,M
        DO j=2,keep+1
        IF(pick1(i)<=odds(j) .AND. pick1(i)>odds(j-1)) THEN
                ma(i)=j-1
        END IF
        IF(pick2(i)<=odds(j) .AND. pick2(i)>odds(j-1)) THEN
                pa(i)=j-1
        END IF
        END DO
        END DO

        !_______________________________________________________
        ! Performs mating using single point crossover

        i=0

        DO i=1,CEILING(REAL(keep)/2.)
        ix(i)=2*i-1
        END DO

        !Allocate temporary array (block 2) (Subroutine requires a real argument)
        ALLOCATE(temp_arr_1(M),STAT=status)
        IF(status/=0) THEN
                WRITE(*,*) 'Error allocating the arrays of allocation block 2 of the main program.'
                STOP
        END IF

        CALL RANDOM_NUMBER(temp_arr_1)

        xp=CEILING(temp_arr_1*REAL(npar))

        DEALLOCATE(temp_arr_1)

        CALL RANDOM_NUMBER(r)

        par2=par

        DO i=1,M
        xy=par2(ma(i),xp(i))-par2(pa(i),xp(i))             ! mix from ma & pa
        par2(keep+ix(i),:)=par2(ma(i),:)                   ! first offspring variable
        par2(keep+ix(i)+1,:)=par2(pa(i),:)                 ! second offspring variable
        par2(keep+ix(i),xp(i))=par2(ma(i),xp(i))-r(i)*xy   ! first offspring variable
        par2(keep+ix(i)+1,xp(i))=par2(pa(i),xp(i))+r(i)*xy ! second offspring variable

        ! Perform crossover when last variable not selected
        IF(xp(i)<npar) THEN                                
                DO j=1,xp(i)
                par2(keep+ix(i),j)=par2(keep+ix(i),j)
                par2(keep+ix(i)+1,j)=par2(keep+ix(i)+1,j)
                END DO

                DO j=xp(i)+1,npar
                par2(keep+ix(i),j)=par2(keep+ix(i)+1,j)
                par2(keep+ix(i)+1,j)=par2(keep+ix(i),j)
                END DO

        END IF
        END DO

        par=par2

        !_______________________________________________________
        ! Mutate the population

        ! Allocate arrays (block 3)
        ALLOCATE(temp_arr_1(nmut),mrow(nmut),mcol(nmut),ind2(nmut),STAT=status)

        IF(status/=0) THEN
                WRITE(*,*) 'Error allocating the arrays of allocation block 3 of the main program.'
                STOP
        END IF

        CALL RANDOM_NUMBER(temp_arr_1)
        mrow=CEILING(temp_arr_1*(popsize-1))+1

        call QsortCInt(mrow)

        CALL RANDOM_NUMBER(temp_arr_1)
        mcol=CEILING(temp_arr_1*npar)

        CALL RANDOM_NUMBER(temp_arr_1)
        temp_arr_1=(hi-lo)*temp_arr_1+lo ! Normalizes values between hi & lo


        DO i=1,nmut
        par(mrow(i),mcol(i))=temp_arr_1(i)
        END DO

        DEALLOCATE(mrow,mcol,temp_arr_1,ind2)

        IF(MINVAL(cost)/=cost(1)) THEN
                bad_sort=bad_sort+1
                IF(bad_sort<=1) THEN
                        WRITE(10,*)cost
                END IF
        END IF

        END DO

        !_______________________________________________________

        DEALLOCATE(par,par2,cost,ind,odds,pick1,pick2,ma,pa,r,xp,ix)
        CLOSE(10)

        WRITE(*,"(I4,' iterations were required to obtain ',I4,' consecutive &
                &values of ',F12.5)") iga,same,minc

END PROGRAM my_cga
