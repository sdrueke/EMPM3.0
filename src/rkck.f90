! -------------------------------------------------------------------------------------------------
! This subroutine takes the Cash-Karp Runge-Kutta step 
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

        SUBROUTINE rkck(y,dydx,n,h,yout,yerr,derivs)

        USE array

        IMPLICIT NONE

        INTEGER       :: i, n
        REAL*8        :: dydx(n), h, y(n), yerr(n), yout(n)

        EXTERNAL derivs

        REAL*8        :: ak2(nmax), ak3(nmax), ak4(nmax), ak5(nmax), ak6(nmax),       &
                      &  ytemp(nmax), A2, A3, A4, A5, A6, B21, B31, B32, B41, B42, B43, B51,  &
                      &  B52, B53, B54, B61, B62, B63, B64, B65, C1, C3, C4, C6, DC1, DC3, DC4, &
                      &  DC5, DC6
        PARAMETER ( A2=0.2,A3=0.3,A4=0.6,A5=1.0,A6=0.875,B21=0.2,B31=3.0/40.0,&
                  & B32=9.0/40.0,B41=0.3,B42=-0.9,B43=1.2,B51=-11.0/54.0,B52=2.5,&
                  & B53=-70.0/27.0,B54=35.0/27.0,B61=1631.0/55296.0,B62=175.0/512.0,&
                  & B63=575.0/13824.0,B64=44275.0/110592.0,B65=253.0/4096.0,C1=37.0/378.0,&
                  & C3=250.0/621.0,C4=125.0/594.0,C6=512.0/1771.0,DC1=C1-2825.0/27648.0,&
                  & DC3=C3-18575.0/48384.0,DC4=C4-13525.0/55296.0,DC5=-277.0/14336.0,&
                  & DC6=C6-0.25 )

        DO i=1,n
           ytemp(i) = y(i)+B21*h*dydx(i)
        END DO

        CALL derivs(ytemp,ak2) !x+A2*h,

        DO i=1,n
           ytemp(i) = y(i)+h*(B31*dydx(i)+B32*ak2(i))
        END DO

        CALL derivs(ytemp,ak3) !x+A3*h,

        DO i=1,n
           ytemp(i) = y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
        END DO

        CALL derivs(ytemp,ak4) !x+A4*h,

        DO i=1,n
           ytemp(i) = y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
        END DO

        CALL derivs(ytemp,ak5) !x+A5*h,

        DO i=1,n
           ytemp(i) = y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
        END DO

        CALL derivs(ytemp,ak6) !x+A6*h,

        DO i=1,n
           yout(i)  = y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
        END DO

        DO i=1,n
           yerr(i)  = h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
        END DO
      RETURN

      END

! -- (C) Copr. 1986-92 Numerical Recipes Software 5"#@-130Rk#3#)KB.
