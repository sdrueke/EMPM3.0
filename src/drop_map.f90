! -------------------------------------------------------------------------------------------------
! This subroutine performs the triplet mapping for droplet inside the linear eddy domain. 
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

	SUBROUTINE drop_map(n,i1,m1,ncell,flag,index,ibound1,ibound2)

        IMPLICIT NONE

        INTEGER :: flag
        INTEGER :: i1, j1, j2, m1, n, nn
        INTEGER :: ncell
        INTEGER :: ibound1, ibound2
        INTEGER, DIMENSION(ibound1:ibound2) :: index

! -- running index
        INTEGER :: j
	
! -- if this is not the first time enter this subroutine, then it is not necessary to 
! -- set the I.C.s and B.C.s again

! -- set I.C.s and B.C.s
	IF (flag .EQ. 0) THEN
	   DO j=1,ncell
	      index(j)       = j
              index(j+ncell) = j
	   END DO
        END IF

! -- Be careful and do not be confused here:
! -- n         amount of eddy grid points
! -- i1-m1+1   first point of triplet mapping (because that index(j)=j, so i1 in LEM is equal 
!              i1-m1+1 in DROP_MAP.f here. 
! -- nn        last point of mapping, is i1-m1+n. 
! -- There should not be a difference, since m1=1. So i1-m1+1=i1
! -- ALso think about if there is difference of X(i) in LEM and DROP)

	CALL triplet1(n,i1-m1+1,index,ibound1,ibound2)
	
	nn = i1-m1+n
!
! -- check if index(j) is negative
	DO j=i1-m1+1,i1-m1+n
           IF (index(j) .LE. 0) THEN
	      WRITE(*,*) n, i1
              WRITE(*,*) 'line 128 m', j, 'index(m)', index(j)
           END IF
	END DO
!
! -- set the first point (j1) 
	j1 = i1-m1+1
	j2 = nn
!
! -- return to index to the initial value
  	DO j=j1,j2
	   IF (j .GT. ncell) THEN
     	      index(j-ncell) = index(j)
	   ELSE 
     	      index(j+ncell) = index(j)
	   END IF
	END DO

        RETURN

        END


! -------------------------------------------------------------------------------------------------
! This subroutine performs the triplet mapping
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

        SUBROUTINE triplet1(kk,k1,A,ibound1,ibound2)

        IMPLICIT NONE

        INTEGER :: b1, k, kk, k1
        INTEGER :: m, ma, mb, mb1, mc, mw
	INTEGER :: ibound1, ibound2
        INTEGER, DIMENSION(ibound1:ibound2) :: A, E

! -- kk       grid number for the eddy, due to the triplet mapping
! -- m        1/3 eddy grid points
!
! -- the logic of triplet mapping can be described as following example
! -- initial scalar: 1,2,3,4,5,6,7,8,9
! -- after mapping:  1,4,7,8,5,2,3,6,9

        m  = kk/3
        ma = 3*m-1
        mb = 3*m-2
        mc = 3*m-3
!
! -- step 1
        DO k=k1,k1+m-1
           E(k) = A(3*m-ma+k1-1)
           ma   = ma-3
        END DO
!
! -- step 2 and step 2.1 - 2nd segment and block inversion
        DO k=m+k1,2*m+k1-1
           E(k) = A(3*m-mb+k1-1)
           mb   = mb-3
        END DO
!
! -- block inversion of the 2nd segment
        mb1 = INT((m-1)/2)+(m+k1)
        mw  = 2*m+k1-1
        DO k=m+k1,mb1
           b1    = E(k)
           E(k)  = E(mw)
        
           E(mw) = b1
           mw    = mw-1
        END DO
!
! -- step3 - 3rd segment
        DO k=2*m+k1,3*m+k1-1
           E(k) = A(3*m-mc+k1-1)
           mc   = mc-3
        END DO

        DO k=k1,kk+k1-1
           A(k) = E(k)
        END DO

        RETURN

        END
