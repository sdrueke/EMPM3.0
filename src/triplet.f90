! -------------------------------------------------------------------------------------------------
! This subroutine performs the triplet mapping for the scalar A
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

	SUBROUTINE triplet(kk,k1,A,igrid)

        IMPLICIT NONE

        INTEGER   :: k, k1, kk, m, ma, mb, mb1, mc, mw
        INTEGER*4 :: igrid
        REAL*8, DIMENSION(0:igrid,2) :: A, E
	REAL*8    :: B1, B2

! -- kk       grid number for the eddy, due to the triplet mapping
! -- m        1/3 eddy grid points
!
! -- the logic of triplet mapping can be described as following example
! -- initial scalar: 1,2,3,4,5,6,7,8,9
! -- after mapping:  1,4,7,8,5,2,3,6,9

        m  = kk/3
        ma = 3*m-1
        mb = 3*m-2
        mc = 3*m-3
!
! -- step 1
        DO k=k1,k1+m-1
           E(k,1) = A(3*m-ma+k1-1,1)
           E(k,2) = A(3*m-ma+k1-1,2)
           ma     = ma-3
        END DO
!
! -- step 2 and step 2.1 - 2nd segment and block inversion
        DO k=m+k1,2*m+k1-1
           E(k,1) = A(3*m-mb+k1-1,1)
	   E(k,2) = A(3*m-mb+k1-1,2)
           mb     = mb-3
        END DO
!
! -- block inversion of the 2nd segment
        mb1 = INT((m-1)/2)+(m+k1)
        mw  = 2*m+k1-1
        DO k=m+k1,mb1
           B1      = E(k,1)
	   B2      = E(k,2)
           E(k,1)  = E(mw,1)
	   E(k,2)  = E(mw,2)
           E(mw,1) = B1
	   E(mw,2) = B2
           mw      = mw-1
        END DO
!
! -- step 3 - 3rd segment 
        DO k=2*m+k1,3*m+k1-1
           E(k,1) = A(3*m-mc+k1-1,1)
	   E(k,2) = A(3*m-mc+k1-1,2)
           mc     = mc-3
        END DO

        DO k=k1,kk+k1-1
           A(k,1) = E(k,1)
	   A(k,2) = E(k,2)
        END DO

        RETURN

        END
