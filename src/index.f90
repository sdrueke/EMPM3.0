! -------------------------------------------------------------------------------------------------
! This function finds the index of the element of table te that is just less than or equal to v
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
! July 24th, 2017
! -- make minor change to enable the use of in-cloud vertical velocity profile which provide only
! -- in-cloud data (starting at cloud base)
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        INTEGER FUNCTION INDEXR( v, ne, te, LF )

! -- glossary
! -- v    variable
! -- ne   number of elements in table te
! -- te   table of elements arranged in ascending or descending order 
! -- LF   logical flag

        IMPLICIT NONE

        INTEGER :: i, j, nd, ne

	REAL*8  :: te(ne),  v
	LOGICAL :: LF

! -- order test
        IF (te(1) .GT. te(2)) GOTO 7

! -- extreme tests
        IF (v .GE. te(1)) GOTO 1
        LF = .TRUE.
        i  = 1
        GOTO 14

   1    IF (v .LT. te(ne)) GOTO 2
        LF = .TRUE.
        i  = ne
        GOTO 14

! -- initializations
   2    LF = .FALSE.
        nd = 1

   3    nd = nd+nd
        IF (nd .GE. ne) GOTO 4
        GOTO 3

   4    nd = nd/2
        i  = nd

! -- bisection loop
   5    nd = nd/2
        IF (nd .LE. 0) GOTO 6
        j = MIN0(ne,i)
        IF (v .GT. te(j)) i = i+nd
        IF (v .LT. te(j)) i = i-nd
        GOTO 5

   6    IF (i .GE. ne) GOTO 14
        IF (v .GT. te(i+1)) i = i+1
        IF (v .LT. te(i))   i = i-1
        GOTO 14

! -- extreme tests
   7    IF (v .GE. te(ne)) GOTO 8
        LF = .TRUE.
        i  = ne
        GOTO 14

   8    IF (v .LT. te(1)) GOTO 9
        LF = .TRUE.
        i  = 1
        GOTO 14

! -- initialization
   9    LF = .FALSE.
        nd = 1

  10    nd = nd+nd
        IF (nd .GE. ne) GOTO 11
        GOTO 10

  11    nd = nd/2
        i  = nd

! -- bisection loop
  12    nd = nd/2
        IF (nd .LE. 0) GOTO 13
        j = MIN0( ne, i )
        IF (v .GT. te(j)) i = i-nd
        IF (v .LT. te(j)) i = i+nd
        GOTO 12

  13    IF (i .GE. ne) GOTO 14
        IF (v .GT. te(i)) i = i-1
        IF (v .LT. te(i)) i = i+1
        i = MAX0(1,i-1)

  14    INDEXR = MIN0(ne-1,i)

        RETURN

        END
