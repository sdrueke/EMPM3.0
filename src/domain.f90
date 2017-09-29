! -------------------------------------------------------------------------------------------------
! This subroutine finds start and end point of the entrained parcel
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

	SUBROUTINE domain(psigma,ngrid,spoint,fpoint,iseed6,n_max,n_blob,n_blob_final,xn)

! -- to decide location to put this parcel into the linear eddy domain.

! -- psigma    portion of the entrainment parcel
! -- xn        grid points of (1-parcel)
! -- spoint    start point of this parcel
! -- fpoint    final point of this parcel

        INTEGER   :: n_max, ngrid
        INTEGER   :: spoint(n_max), fpoint(n_max), n_s(n_max), n_f(n_max)
        INTEGER   :: n_blob, n_blob_final
        INTEGER   :: i, k, xn
	INTEGER*4 :: iseed6
	REAL*8    :: psigma, rand_no
        REAL*8    :: RAND1

! -- calculate how many grid points (1-parcel) has
        WRITE(7,*) ''
	WRITE(7,*) '---------------------------------------------------------------------------'
        write(7,*) 'Starting with domain'

        xn = INT((1.0-psigma*n_blob)*ngrid)

        WRITE(7,*) 'psigma = ', psigma
        WRITE(7,*) 'ngrid  = ', ngrid
        WRITE(7,*) 'xn     = ', xn

! -- randomly choosing the spoint and fpoint shown by grid point relative to this parcel.
	
	rand_no = RAND1(iseed6)
        i = 1
        n_s(i) = INT(rand_no*xn)+1
        n_f(i) = n_s(i)+psigma*ngrid-1
        DO i=2,n_blob
           n_s(i) = n_f(i-1)+INT(xn/n_blob) 
           n_f(i) = n_s(i)+psigma*ngrid-1
        END DO
        
! -- check if n_s or n_f are greater than ngrid, if it is, it needs to re-establish period
! -- condition
        k = 0
        DO i=1,n_blob
           k = k+1
           if( (n_s(i) .LE. ngrid) .AND. (n_f(i) .GT. ngrid) ) THEN
                spoint(k)   = n_s(i)
                fpoint(k)   = ngrid
                spoint(k+1) = 1
                fpoint(k+1) = n_f(i)-ngrid
                k = k + 1 ! increase k by 1
           ELSE IF ((n_s(i) .GT. ngrid) .AND. (n_f(i) .GT. ngrid) ) THEN
                spoint(k) = n_s(i)-ngrid
                fpoint(k) = n_f(i)-ngrid
           ELSE
                spoint(k) = n_s(i)
                fpoint(k) = n_f(i)
           END IF
        END DO

        n_blob_final = k 

        WRITE(7,*) '          k       spoint      fpoint'
	WRITE(7,*) '---------------------------------------------------------------------------'
        DO k=1,n_blob_final
           WRITE(7,*) k, spoint(k), fpoint(k)
        END DO 
    
        WRITE(7,*) 'n_blob_final = ', n_blob_final

        WRITE(7,*) 'Leaving domain'
	WRITE(7,*) '---------------------------------------------------------------------------'

        RETURN
        END
