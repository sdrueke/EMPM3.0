! -------------------------------------------------------------------------------------------------
! This subroutine calculates the droplet growth (droplet growth model - DGM) by Dr. Phil Austin
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

        SUBROUTINE DGM(radius,qv,temp,SuS,press,ver_vel,height,ql,weight,i_grid_vol,t_use,ite,ncat)
     
! --    array drop_radius:  1--ndmax contain the radii of each droplet category
!     
! --    the  other variables:
! --    ndmax+1 = wv water vapor (kg/kg)
! --    ndmax+2 = temp (K)
! --    ndmax+3 = s (fraction)
! --    ndmax+4 = p pressure (pa)
! --    ndmax+5 = w constant vertical velocity (m/s)
! --    ndmax+6 = z 
! --    ndmax+7 = wl liquid water (kg/kg)
        USE const
        USE array

        IMPLICIT NONE

        EXTERNAL fcnkb,rkqs    
 
        INTEGER :: ite
! -- aerosol category (ncat)
        INTEGER :: nbad, ncat, nok

        REAL*8  :: i_grid_vol
        REAL*8  :: height, press, SuS, ver_vel
        REAL*8  :: epserr, es, ew, h1, hmin
        REAL*8  :: ql, qs, qv 
        REAL*8  :: radius
        REAL*8  :: SS, t1, t2
        REAL*8  :: temp, t_use 
        REAL*8  :: weight

	REAL*8, DIMENSION(nmax) :: drop_radius

 
! --    SuS     = Supersaturation
! --    press   = pressure (pa)
! --    ver_vel = w, constant velocity
! --    height  = height
! --    ql      = liquid water mixing ratio (kg/kg)

	ndmax    = 1
	dmaxa    = ncat

        drop_radius(1) = radius
        grid_scale     = i_grid_vol
        solute_mass    = weight

        drop_radius(2) = qv
        drop_radius(3) = temp
        drop_radius(5) = press
        drop_radius(6) = ver_vel
        drop_radius(7) = height

	IF (ite .LE. 1) THEN
           drop_radius(8) = pi43*rho_w*drop_radius(1)**3*grid_scale
        ELSE
           drop_radius(8) = ql
        END IF
 
        es       = ew(drop_radius(3))*100.0
        qs       = eps*es/(drop_radius(5)-es)
        SS       = drop_radius(2)/qs-1.0

        drop_radius(4) = SS

        t1       = 0.0
        t2       = 0.0
        epserr   = 1.0e-4
        hmin     = 0.0

! -- add time steps for odeint
        h1       = t_use*1.0
        t1       = ite*t_use
        t2       = t1+h1

        CALL odeint(drop_radius,8,t1,t2,epserr,h1,hmin,nok,nbad,fcnkb,rkqs)

        radius   = drop_radius(1)

        qv       = drop_radius(2)
        temp     = drop_radius(3)
        SuS      = drop_radius(4)
        press    = drop_radius(5)
        ver_vel  = drop_radius(6)
        height   = drop_radius(7)
        ql       = drop_radius(8)

        RETURN

        END
