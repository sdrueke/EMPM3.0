! -------------------------------------------------------------------------------------------------
! This subroutine takes the fifth-order Runge-Kutta step with monitoring of local truncation error
! to ensure accuracy and adjust stepsize
! -------------------------------------------------------------------------------------------------
        SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)

        USE array

        IMPLICIT NONE

        INTEGER       :: i, n
        REAL*8        :: eps, errcon, errmax, h, hdid, hnext, htry, pgrow, pshrnk, safety, x, xnew
        REAL*8        :: dydx(n), y(n), yerr(nmax), yscal(n), ytemp(nmax)
        PARAMETER (safety=0.9,pgrow=-0.2,pshrnk=-0.25,errcon=1.89e-4)

        EXTERNAL derivs

        h      = htry

1       CALL rkck(y,dydx,n,h,ytemp,yerr,derivs)

        errmax = 0.0
        DO i=1,n
           errmax = MAX(errmax,ABS(yerr(i)/yscal(i)))
        END DO

        errmax = errmax/eps

        IF (errmax .GT. 1.0) THEN
           h = safety*h*(errmax**pshrnk)
           IF (h .LT. 0.1*h) THEN
              h = 0.1*h
           END IF
           xnew = x+h
           IF (xnew .EQ. x) THEN
              WRITE(7,*) 'stepsize underflow in rkqs'
              STOP
           END IF
           GOTO 1
        ELSE
           IF(errmax .GT. errcon) THEN
              hnext = safety*h*(errmax**pgrow)
           ELSE
              hnext = 5.0*h
           END IF
           hdid = h
           x    = x+h
           DO i=1,n
              y(i) = ytemp(i)
	      IF ((i .LE. ndmax) .AND. (y(i) .LE. 0.0)) y(i) = 1.d-8
           END DO
           RETURN
        END IF

        END

! -- (C) Copr. 1986-92 Numerical Recipes Software 5"#@-130Rk#3#)KB.
