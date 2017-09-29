! -------------------------------------------------------------------------------------------------
! This subroutine calculates the accumulated concerntration of the initial droplet distribution
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! April 22nd, 2016
! -- bug fix: entrained ccn read from file works properly now
! -- simplified the droplet number calculation in the domain
! -- removed the redistribution of the droplets (original distribution of the input file is used)
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

	SUBROUTINE accu(volume,dis,con_drop,ndrop,ne,ne1)

        IMPLICIT NONE

        INTEGER :: i
	INTEGER :: ndrop, ne, ne1

        REAL*8  :: volume

        REAL*8  :: dis(ne)
	REAL*8  :: con_drop(ne1)

        WRITE(7,*) ''
	WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) 'Starting with accu'

        DO i=1,ne
           dis(i) = 0.0
        END DO

        DO i=1,ne
           dis(i)     = NINT(con_drop(i)*volume)
        END DO

	WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) ' i   distribution    droplet concentration'
        WRITE(7,*) '      (#/domain)      (#/m^3)' 

        DO i=1,ne
           WRITE(7,100) i, dis(i), con_drop(i)
        END DO
100 format(i3,3x,e12.5,4x,e12.5)

        ndrop = SUM(INT(dis))
!        WRITE(71,*)'number of drop is', ndrop

        WRITE(7,*) 'Leaving accu'
	WRITE(7,*) '---------------------------------------------------------------------------'

	RETURN

	END
