! -------------------------------------------------------------------------------------------------
! This function interpolates between X1 and X2 using a functional form that depends upon the 
! specified mode
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

        REAL*8 FUNCTION FINTRP ( MODE, X, X1, F1, X2, F2 )

	REAL*8 :: F, F1, F2, X, X1, X2

        INTEGER :: MODE
! -- 1     linear       -- A * X + B
! -- 2     exponential  -- A * EXP( B * X )
! -- 3     power law    -- A * X**B

        GOTO ( 10, 20, 30 ), MODE

! -- linear interpolation
  10    F = F1+(F1-F2)*(X-X1)/(X1-X2)
        GOTO 40

! -- exponential interpolation
  20    F = F1*(F1/F2)**((X-X1)/(X1-X2))
        GOTO 40

! -- power law interpolation
  30    IF (X1 .LE. 0.0) GOTO 20
        F = F1*(X/X1)**(DLOG(F1/F2)/DLOG(X1/X2))

  40    FINTRP = F

        RETURN

        END
