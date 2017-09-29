! -------------------------------------------------------------------------------------------------
! This subroutine calculates the pdf of the scalar. There are 100 bins.
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

        SUBROUTINE prob(m1,m4,pNmin,pNmax,A,PDF,SCALAR,scalar_pdf_nbin,ipe,mip,igrid)

! -- PDF      power density function
! -- SCALAR   scalar values divided by 100 intervals
! -- nbin     bin numbers for PDF
! -- PDFAREA  the total area of PDF
! -- j        locate the A(k) into approciate bin #
! -- pNmin    minimum of scalar value
! -- pNmax    maximum of scalar value
! -- ipe      every period of one realization, if periods=25 iteration
!             then they are 25, 50, 75, 100...... iteration for periods 1,
!             2, 3, 4, ......

        IMPLICIT NONE

        INTEGER   :: ipe, j, k, m1, m4, mip, scalar_pdf_nbin, nbin
        INTEGER*4 :: igrid
        REAL*8    :: ddx, pNlength, pNmax, pNmin, PDFAREA
        REAL*8    :: A(0:igrid), PDF(mip,100), SCALAR(100)

        nbin = scalar_pdf_nbin

! -- set initial PDF values zero
        DO k=1,nbin
           PDF(ipe,k) = 0.0
        END DO

        pNlength = pNmax-pNmin
        PDFAREA  = 0.0

        DO k=m1,m4
           j = INT(nbin*((A(k)-pNmin)/pNlength))+1
           IF (j .GE. 101 .OR. j .LE. 0) THEN
              WRITE(7,*) 'pNmin, pNmax =', pNmin, ' ',  pNmax
              WRITE(7,*) 'j = ', j, ' a(k)= ', a(k)
              j = 100
           END IF
! -- accumulate the count # at each interval of j
           PDF(ipe,j) = PDF(ipe,j)+1.0
        END DO

! -- calculate the total PDFAREA
        ddx = (pNlength)/FLOAT(nbin)

        DO k=1,nbin
           PDFAREA = PDFAREA+PDF(ipe,k)*ddx
        END DO

! -- normlaized the PDF and calculate the real scalar value at each bin#
        DO k=1,nbin
           PDF(ipe,k) = PDF(ipe,k)/PDFAREA
           SCALAR(k)  = (k-1.0)*ddx+pNmin+ddx/2.0
        END DO

        RETURN

        END



