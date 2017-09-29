! -------------------------------------------------------------------------------------------------
! This subroutine calculates the droplet size distribution using a gamma function
! -------------------------------------------------------------------------------------------------
!
! Date:        March 2016
!
!
! Updates:
! --------
!
! May 18th, 2016
! -- fixing a bug which led to an overestimation of liquid water mixing ratio
!
! March 7th, 2016
! -- code was created
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------
        SUBROUTINE size_dis(N_i,qc_i,rho,radius_drop,con_drop,mass_drop,volume,m_i,ne)

        USE const

        IMPLICIT NONE

        INTEGER :: n
        INTEGER :: ne
        REAL*8  :: N_i, volume, rho
        REAL*8  :: qc_i, m_i

        REAL*8  :: mu, lambda
        REAL*8  :: const1, const2, iconst2

        REAL*8  ::  N0
        REAL*8, DIMENSION(ne) :: diameter_drop, delta_diameter
        REAL*8, DIMENSION(ne) :: radius_drop, con_drop, mass_drop
        REAL*8  :: min_diameter = 1.0E-6

        REAL*8  :: gamma1, gamma4, igamma1

        WRITE(7,*) ''
        WRITE(7,*) '---------------------------------------------------------------------------'
        WRITE(7,*) 'size_dis.f90: (calculation of size distribution using gamma function)'

        const1  = pi*rho_w/6.0
        const2  = 3.0
        iconst2 = 1.0/const2

! -- Calculate shape paremeter of the gamma function
! -- Morrison and Gettelman,, 2008
        mu = 0.0005714*N_i+0.2714
        mu = 1./(mu**2.)-1.
        mu = MAX(mu,2.)
        mu = MIN(mu,10.)      

! -- Thompson et al, 2008
        !mu = MIN(15.,1000.E6/N_i+2)

        WRITE(7,*) 'mu     = ', mu

        gamma1 = GAMMA(mu+1)
        gamma4 = GAMMA(mu+const2+1)

        igamma1 = 1.0/gamma1

! -- Calculate distribution's slope
        lambda = 1.0D-6*(N_i*rho*const1*gamma4*igamma1/qc_i)**iconst2

        WRITE(7,*) 'lambda = ', lambda

        N0 = 1.0D-18*N_i*lambda**(mu+1)*igamma1

! -- Create bins of cloud water 
        diameter_drop(1)  = min_diameter*1.0d0
        delta_diameter(1) = min_diameter*1.0d0
        DO n = 2,ne
           diameter_drop(n)  = diameter_drop(n-1)+1.0D-6
           delta_diameter(n) = (diameter_drop(n)-diameter_drop(n-1))
        END DO

! -- Calculate size distribution
        DO n = ne,1,-1
           diameter_drop(n) = diameter_drop(n)*1.0D6
           con_drop(n)      = N0*diameter_drop(n)**mu* EXP(-lambda*diameter_drop(n))*delta_diameter(n)
           con_drop(n)      = 1.0D24*con_drop(n)*volume
           radius_drop(n)   = (diameter_drop(n)*1.0D-6)/2
        END DO

        DO n=1,ne
           IF (NINT(con_drop(n)) .NE. 0.0) THEN
              mass_drop(n) = m_i 
           END IF
        END DO
 
        END SUBROUTINE
