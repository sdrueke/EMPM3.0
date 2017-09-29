! -------------------------------------------------------------------------------------------------
! This subroutine calculates the effective radius
! -------------------------------------------------------------------------------------------------
!
! Date:        August 2016
!
!
! Updates:
! --------
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        SUBROUTINE eff_radius(n_drop,nbin,r2_mean,r3_mean,pdf,r_eff)

        IMPLICIT NONE

        INTEGER :: i, nbin, n_drop
        REAL*8  :: numerator, denominator
        REAL*8  :: r2_mean(nbin), r3_mean(nbin), pdf(nbin)
        REAL*8  :: r_eff

        numerator   = 0.0
        denominator = 0.0

        DO i=1,nbin
           pdf(i)      = pdf(i)*n_drop

           numerator   = numerator+(pdf(i)*r3_mean(i))
           denominator = denominator+(pdf(i)*r2_mean(i))
        END DO

        r_eff = numerator/denominator

        RETURN
        END

