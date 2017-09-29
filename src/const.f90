! -------------------------------------------------------------------------------------------------
! This module sets constants 
! -------------------------------------------------------------------------------------------------
        MODULE const

        REAL*8, PARAMETER :: pi   = 3.1415926535897931d+0
        REAL*8, PARAMETER :: pi4  = 12.566370614359172d+0
        REAL*8, PARAMETER :: pi43 = 4.1887902047863905d+0
! gravitational constant (g), dynamic viscosity (dyn_vis)
	REAL*8, PARAMETER :: g = 9.81, dyn_vis=1.488e-5
! specific heat of dry air at constant pressure (cp) 
! specific heat of dry air at constant pressure (cv) 
! specific heat of water vapour at constant pressure (cpv)
        REAL*8, PARAMETER :: cp = 1005.0, cv = 718.0
        REAL*8, PARAMETER :: cpv = 1875.0
! specific heat of water (cw)
        REAL*8, PARAMETER :: cw = 4190.0
! constant latent heat (Lv_c)
        REAL*8, PARAMETER :: Lv_c = 2.5e6
! individual gas constant of water vapour (Rv)
! individual gas constant of dry air (Rd)
        REAL*8, PARAMETER :: Rv = 461.5, Rd = 287.0
! universal gas constant (R_uni)
        REAL*8, PARAMETER :: R_uni = 8.1344598
! molecular weight of air (Md) and of water (Mw)
        REAL*8, PARAMETER :: Md = 28.96535!E-3
        REAL*8, PARAMETER :: Mw = 18.01528

        REAL*8, PARAMETER :: eps = 0.622, epstv = 0.608 
! and water density (rho_w)
        REAL*8, PARAMETER :: rho_w = 1000.0

! temperature dependent latent heat (Lv)
        REAL*8 :: Lv
! diffusivity
        REAL*8 :: Ktemp, D
! specific heat of air, water vapor and liquid water
        REAL*8 :: cpm

        SAVE

        END MODULE const
