! -------------------------------------------------------------------------------------------------
! This subroutine calculates the maximum of droplet radius
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        SUBROUTINE rmax(datamax,data,ndrop,index_rmax)

        IMPLICIT NONE

        INTEGER :: index_rmax, j, ndrop
        REAL*8  :: datamax, datadiff
        REAL*8  :: data(ndrop)

! -- calculate maximum
        DO j=1,ndrop
           datadiff = data(j)*1.e6-datamax
           IF (datadiff .GT. 0.0) THEN
              datamax    = data(j)*1.e6
	      index_rmax = j
           END IF
        END DO

        RETURN

        END
