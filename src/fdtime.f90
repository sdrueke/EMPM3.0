! -------------------------------------------------------------------------------------------------
! This subroutine calculates the time interval, t_prime, between entrainment events
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! August 3rd, 2016
! -- combine two parts of EMPM (both calculated entrainment time step) and simply code
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        SUBROUTINE fdtime(psigma,n_blob,ent_rate,iseed3,w,ire,dt_entm)

        IMPLICIT NONE

        INTEGER   :: ire, n_blob
        INTEGER*4 :: iseed3
        REAL*8    :: c, dem, delta_z, ent_rate
        REAL*8    :: RAND1
	REAL*8    :: dt_entm, psigma, w 

        WRITE(7,*) ''
	WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) 'Starting with fdtime'

! -- determine entrainment frequency (time interval between two entrainment events)
        delta_z  = (n_blob/ent_rate)*(psigma/(1.0-psigma))
        dt_entm  = delta_z/ABS(w)

        IF (ire .EQ. 1) THEN

! -- c is a random # between 0,1; DLOG(1./1-c) will have a mean time 1.
! -- and an exponetial distribution, Poisson distribution
           c       = RAND1(iseed3)
           dem     = -DLOG(1.-c)
! -- this is for ave. time of entm event
           dt_entm = dt_entm*dem

        END IF

        WRITE(7,*) 'dt_entm = ', dt_entm
        write(7,*) 'Leaving fdtime'
	WRITE(7,*) '---------------------------------------------------------------------------'

        RETURN

        END

