! -------------------------------------------------------------------------------------------------
! This module sets variables for spline interpolation 
! -------------------------------------------------------------------------------------------------
        MODULE intspline

! dimension for splinearray 
        INTEGER, PARAMETER :: ncats = 49
        INTEGER, PARAMETER :: nvecs = 3
        INTEGER, PARAMETER :: npoints = 11

! array for spline interpolation (splinearray)
        REAL*8 :: splinearray(ncats,nvecs,npoints)

        SAVE

        END MODULE
