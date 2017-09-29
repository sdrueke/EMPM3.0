! -------------------------------------------------------------------------------------------------
! This subroutine randomly assigns the droplet position
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

	SUBROUTINE assigd(ndrop,ngrid,jcell_pop,dx,m1,x,BL,igrid,idimen,iseed4)

        IMPLICIT NONE

        INTEGER   :: i, j, m1, idimen
        INTEGER   :: ndrop, ngrid
        INTEGER   :: nis1, nis2, nis3, nis4, ngt4
	INTEGER   :: jcell_pop(0:10,igrid)
	INTEGER*4 :: igrid, iseed4
        REAL*8    :: BL
        REAL*8    :: RAND1
	REAL*8    :: x(idimen)
	REAL*8    :: dx, rand_no

        WRITE(7,*) 'Starting with assigd'

! -- initialize jcell_pop(0,j) which is the index used in drop_map.f90
	DO j=1,ngrid
	   jcell_pop(0,j) = 0
	END DO

        DO i=1,ndrop
	   rand_no = RAND1(iseed4)

	   x(i)    = rand_no*BL
	   j       = INT(x(i)/dx)+M1
	   IF (j .GT. ngrid) j = j-ngrid

	   jcell_pop(0,j-m1+1)                   = jcell_pop(0,j-m1+1)+1
	   jcell_pop(jcell_pop(0,j-m1+1),j-m1+1) = i
	END DO

! -- write out how many cells have one, two, three etc droplets (just for fun)
	ngt4 = 0
        nis4 = 0
	nis3 = 0
        nis2 = 0
        nis1 = 0
	
        DO j = 1, ngrid
	   IF (jcell_pop(0,j-m1+1) .GT. 4) THEN
	      ngt4 = ngt4 + 1
           ELSE IF (jcell_pop(0,j-m1+1) .GT. 3) THEN
	      nis4 = nis4 + 1
	   ELSE IF (jcell_pop(0,j-m1+1) .GT. 2) THEN
	      nis3 = nis3 + 1
	   ELSE IF (jcell_pop(0,j-m1+1) .GT. 1) THEN
	      nis2 = nis2 + 1
	   ELSE IF (jcell_pop(0,j-m1+1) .GT. 0) THEN
	      nis1 = nis1 + 1
	   END IF      
        END DO

        WRITE(7,*) ngt4, 'cells have >4 droplets'
	WRITE(7,*) nis4, 'cells have 4 droplets'
        WRITE(7,*) nis3, 'cells have 3 droplets'
        WRITE(7,*) nis2, 'cells have 2 droplets'
        WRITE(7,*) nis1, 'cells have 1 droplet'

        WRITE(7,*) 'Leaving assigd'
	WRITE(7,*) '---------------------------------------------------------------------------'

	RETURN
	END

! -------------------------------------------------------------------------------------------------
! This subroutine randomly assigns the entrained ccn locations within m2 and m3, which are the 
! first and last grid cells for the newly entrained blob
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! July 19th, 2016
! -- change variable type for RAND1 to REAL*8 - fixed bug related to determination of indices of 
! -- droplet
!
! January 25th, 2016
! -- the old FORTRAN77 code was rewritten in FORTRAN90
! -- full control over code was gained by declaring each variable
! -- ATTENTION: The code does not work!
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

	SUBROUTINE assigd2(ndrop_used,nccn,m2,m3,jcell_pop,dx,m1,x,igrid,idimen,iseed5,n_max, &
                   &       n_blob,part,blob_flag,ngrid)

        IMPLICIT NONE

        INTEGER   :: idimen, ndrop_used, nccn, m1, n_max, flag, part, blob_flag, ngrid
        INTEGER   :: n_blob, n_blob_total
        INTEGER, DIMENSION(n_max) :: m2, m3
	INTEGER   :: jcell_pop(0:10,igrid)
        INTEGER   :: i, j, k
        INTEGER   :: ngt3, nis1, nis2, nis3
	INTEGER*4 :: igrid
	INTEGER*4 :: iseed5
        REAL*8    :: RAND1
        REAL*8    :: x(idimen)
	REAL*8    :: dx, rand_no

        IF (blob_flag .EQ. 1) THEN
           n_blob_total = n_blob+1
        ELSE
           n_blob_total = n_blob
        END IF

! -- initialize jcell_pop(0,j) which is the index used in drop_map.f90
        DO k=1,n_blob_total
	   DO j=m2(k),m3(k)
	      jcell_pop(0,j-m1+1) = 0
           END DO
	END DO

        flag = 1
        k = 1
	DO i=ndrop_used+1,ndrop_used+nccn
! -- randomly assign x(i)
! -- find corresponding j_cell for x(i)
	      rand_no = RAND1(iseed5)

              IF (blob_flag .EQ. 1 .AND. k .EQ. n_blob) THEN
                 IF (flag .LE. part) THEN
                    flag = flag+1
                 ELSE IF (flag .GT. part .AND. flag .LT. 10) THEN
                    k    = n_blob+1
                    flag = flag+1
                 ELSE
                    k    = n_blob+1
                    flag = 1
                 END IF
              END IF

	      x(i) = rand_no*(m3(k)-m2(k)+1)*dx+(m2(k)-1)*dx
	      j    = INT(x(i)/dx)+M1
	      IF (j .GT. ngrid) j = j-ngrid

              jcell_pop(0,j-m1+1)                   = jcell_pop(0,j-m1+1)+1
 	      jcell_pop(jcell_pop(0,j-m1+1),j-m1+1) = i-ndrop_used

              IF (k .GE. n_blob) THEN
                 k = 1
              ELSE
                 k = k+1
              END IF
	END DO
	
! -- write out how many cells have one, two, three etc droplets (just for fun)
        ngt3 = 0
        nis3 = 0
        nis2 = 0
        nis1 = 0



        DO k=1,n_blob_total
           DO j = m2(k), m3(k)
              IF (jcell_pop(0,j-m1+1) .GT. 3) THEN
                 ngt3 = ngt3 + 1
              ELSE IF (jcell_pop(0,j-m1+1) .GT. 2) THEN
                 nis3 = nis3 + 1
              ELSE IF (jcell_pop(0,j-m1+1) .GT. 1) THEN
                 nis2 = nis2 + 1
              ELSE IF (jcell_pop(0,j-m1+1) .GT. 0) THEN
                 nis1 = nis1 + 1
              END IF
           END DO
        END DO

	WRITE(7,*) ngt3, 'cells have >3 droplets'
        WRITE(7,*) nis3, 'cells have 3 droplets'
        WRITE(7,*) nis2, 'cells have 2 droplets'
        WRITE(7,*) nis1, 'cells have 1 droplet'

	RETURN

	END

