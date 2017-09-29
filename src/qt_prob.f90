! -------------------------------------------------------------------------------------------------
! This subroutine calculates the probability density function for total water mixing ration for a 
! movie???
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

        SUBROUTINE qt_prob(m1,m4,pNmin,pNmax,A,PDF,SCALAR,igrid,qt_pdf_nbin)

        IMPLICIT NONE

        INTEGER :: j, k, m1, m4, nbin
	INTEGER :: igrid, qt_pdf_nbin
        REAL*8  :: A(igrid)
        REAL*8  :: PDF(qt_pdf_nbin), SCALAR(qt_pdf_nbin)
	REAL*8  :: ddx, pNlength, pNmin, pNmax, PDFAREA

        nbin = qt_pdf_nbin

        DO k=1,nbin
           PDF(k) = 0.0
        END DO

! -- PDFAREA  the total area of PDF
! -- k        locate the A(k) into approciate bin #
! -- pNmin    minimum of scalar value
! -- pNmax    maximum of scalar value

        pNlength = pNmax-pNmin
        PDFAREA  = 0.0

        DO k=m1,m4
           j = INT(nbin*((A(k)-pNmin)/pNlength))+1

           IF(j .GT. nbin .OR. j .LE. 0) THEN
!              WRITE(7,*) 'qt_prob.f90: j, a(k) = ', j, a(k)
              j = nbin
           END IF
! -- accumulate the count # at each interval of j
           PDF(j) = PDF(j)+1.0
        END DO

! -- calculate the total PDFAREA
        ddx = (pNlength)/FLOAT(nbin)

        DO k=1,nbin
           PDFAREA = PDFAREA+PDF(k)*ddx
        END DO

! -- normlaized the PDF and calculate the real scalar value at each bin#
        DO k=1,nbin
           PDF(k)    = PDF(k)/PDFAREA
           SCALAR(k) = (k-1.0)*ddx+pNmin+ddx/2.0
        END DO

        RETURN
        END

