! -------------------------------------------------------------------------------------------------
! This subroutine interpolates the environmental data given to the pressure level pe
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! June 9th, 2016
! -- remove pressure calculation and use the in the droplet growth model calculated pressure 
!    instead
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
       ! SUBROUTINE intpol(time0,w,ni,p1,vi,pe,v)
        SUBROUTINE intpol(ni,p1,vi,pe,v)
        USE const

        IMPLICIT NONE

! -- pe      pressure level to interpolate to
! -- pn      pressure level in hPa
! -- w       vertical velocity of the parcel
! -- vi(k)   array of sounding variable with ni elements
! -- v       (lineraly) interpolated variable at pressure level pe

        INTEGER :: k1, k2, ni
        INTEGER :: INDEXR
      !  REAL*8  :: time0
        REAL*8  :: FINTRP
        REAL*8  :: p1(ni), vi(ni)
	REAL*8  :: pe, pn, v!, w
	LOGICAL :: LF

        pn = pe*0.01

        k1 = INDEXR ( pn, ni, p1, LF )
        k2 = k1 + 1

        v = FINTRP ( 1, pn, p1(k1), vi(k1), p1(k2), vi(k2) )

        RETURN
        END

