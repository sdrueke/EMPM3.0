! -------------------------------------------------------------------------------------------------
! This subroutine is the Runge-Kutta driver with adaptive stepsize control
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

        SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)

        USE array

        IMPLICIT NONE

        INTEGER       :: i, kmax, kmaxx, kount, maxstp
        INTEGER       :: nbad, nok, nstp, nvar  
        REAL*8        :: eps, h1, hmin, x1, x2, TINY
        PARAMETER (maxstp=10000,kmaxx=200,TINY=1.e-30)
        REAL*8        :: dydx(nmax), ystart(nvar)
        REAL*8        :: dxsav, h, hdid, hnext, x, xsav
        REAL*8        :: xp(kmaxx), y(nmax), yp(nmax,kmaxx), yscal(nmax)

        EXTERNAL derivs,rkqs

        dydx(:) = 0.0

        x     = x1
        h     = SIGN(h1,x2-x1)
        nok   = 0
        nbad  = 0
        kount = 0
        kmax  = 0

        DO i=1,nvar
           y(i) = ystart(i)
        END DO

        IF (kmax .GT. 0) xsav=x-2.*dxsav

        DO nstp=1,maxstp
           CALL fcnkb(y,dydx)
           DO i=1,nvar
              yscal(i) = ABS(y(i))+ABS(h*dydx(i))+TINY
           END DO

           IF (kmax .GT. 0) THEN
              IF (ABS(x-xsav) .GT. ABS(dxsav)) THEN
                 IF (kount .LT. kmax-1) THEN
                    kount     = kount+1
                    xp(kount) = x
                    DO i=1,nvar
                       yp(i,kount) = y(i)
                    END DO
                    xsav = x
                 END IF
              END IF
           END IF

           IF ((x+h-x2)*(x+h-x1) .GT. 0.0) h=x2-x

           CALL rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)

           IF (hdid .EQ. h) THEN
              nok  = nok+1
           ELSE
              nbad = nbad+1
           END IF

           IF ((x-x2)*(x2-x1) .GE. 0.0) THEN
              DO i=1,nvar
                 ystart(i)=y(i)
              END DO

              IF (kmax .NE. 0) THEN
                 kount=kount+1
                 xp(kount)=x
                 DO i=1,nvar
                    yp(i,kount)=y(i)
                 END DO
              END IF
              RETURN
           END IF

           IF (ABS(hnext) .LT. hmin) THEN
              WRITE(7,*) 'stepsize smaller than minimum in odeint'
              STOP
           END IF

           h = hnext
        END DO
        WRITE(7,*) 'too many steps in odeint'
        RETURN

        END

! --  (C) Copr. 1986-92 Numerical Recipes Software 5"#@-130Rk#3#)KB.
