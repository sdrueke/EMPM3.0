! -------------------------------------------------------------------------------------------------
! This subroutine carries out the spline interpolation for very small droplet 
! -------------------------------------------------------------------------------------------------
!
! Date:        September 2017
!
! Updates:
! --------
!
! September 21st, 2017
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

      SUBROUTINE psplint(ncat,x,y)
      !include "interp_f90.h"

      USE intspline

      IMPLICIT NONE

      INTEGER :: ncat, klo, khi, k
      REAL*8 h,a,b,y,x

! -- this started out as numerical recipies splint, uses derivitives and points passed in splinearray
! -- through intspline.f90
      
      IF (x .LT. splinearray(ncat,1,1) .OR. x .GT. splinearray(ncat,1,npoints)) THEN
         WRITE(7,*) 'interpolation out of range, for category = ',ncat
         WRITE(7,*) 'asking for: ', x, 'range is: ', splinearray(ncat,1,1)
         WRITE(7,*) 'acutal value is: ',  splinearray(ncat,1,npoints)
         STOP
      END IF

      klo = 1
      khi = npoints

 1    IF (khi-klo .GT. 1) THEN
         k = (khi+klo)/2
         IF (splinearray(ncat,1,k) .GT. x) THEN
            khi = k
         ELSE
            klo = k
         END IF
         GOTO 1
      END IF

      h = splinearray(ncat,1,khi)-splinearray(ncat,1,klo)
!c     if (h.eq.0.) pause 'bad splinearray input.'
      a = (splinearray(ncat,1,khi)-x)/h
      b = (x-splinearray(ncat,1,klo))/h
      y = a*splinearray(ncat,2,klo)+b*splinearray(ncat,2,khi)+ &
          ((a**3-a)*splinearray(ncat,3,klo)+                   &
          (b**3-b)*splinearray(ncat,3,khi))*(h**2)/6.0

      RETURN

      END
