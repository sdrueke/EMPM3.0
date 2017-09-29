! -------------------------------------------------------------------------------------------------
! This subroutine finds the domain mean and rms of qlwater, qw, sl, wdrop, hm, SuS, qv and T
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

        SUBROUTINE mean_rms(m1,m4,scalar,scalar_mean,scalar_rms,ip,mip,igrid)

        IMPLICIT NONE

        INTEGER   :: i, ip, j, m1, m4, mip, n
	INTEGER*4 :: igrid
        REAL*8    :: scalar(1:8,igrid)
	REAL*8    :: scalar_mean(1:8,0:mip), scalar_rms(1:8,0:mip)
        REAL*8    :: scalar_sum(1:8), scalar_prim(1:8)

	scalar_sum(:)  = 0.0
	scalar_prim(:) = 0.0

	DO i=m1,m4
           DO j=1,8
	      scalar_sum(j) = scalar_sum(j)+scalar(j,i)
           END DO
	END DO
	
	n = m4-m1+1
        DO j=1,8
	   scalar_mean(j,ip) = scalar_sum(j)/FLOAT(n)
        END DO

	DO i=m1,m4
           DO j=1,8
	      scalar_prim(j) = scalar_prim(j)+(scalar(j,i)-scalar_mean(j,ip))**2.0
           END DO
	END DO

        DO j=1,8
	   scalar_rms(j,ip) = (scalar_prim(j)/FLOAT(n))**0.5
        END DO

	RETURN 
	
        END
