! -------------------------------------------------------------------------------------------------
! This subroutine calculates the maximum and minimum of input variables
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

	SUBROUTINE max_min(datamax,datamin,data,m1,m4,igrid)

        IMPLICIT NONE

        INTEGER   :: j, m1, m4
	INTEGER*4 :: igrid
	REAL*8    :: data(igrid)
        REAL*8    :: datamax, datamin, datadiff

! -- calculate maximum
        DO j=m1,m4
           datadiff = data(j)-datamax
           IF (datadiff .GT. 0.0) THEN
              datamax = data(j)
	   END IF
	END DO

! -- calculate minimum
	DO j=m1,m4
           datadiff = data(j)-datamin
           IF (datadiff .LT. 0.0) THEN
              datamin = data(j)
           END IF
        END DO

	RETURN

        END
