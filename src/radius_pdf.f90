! -------------------------------------------------------------------------------------------------
! This subroutine calculates the probability density function for droplet radius
! -------------------------------------------------------------------------------------------------
!
! Date:        January 2016
!
!
! Updates:
! --------
!
! April 28th, 2016
! -- combined radius_pdf.f90 and mov_sort.f90 in radius_pdf.f90 (deleted mov_sort.f90)
!
! January 28th, 2016
! -- included subroutine to calculate the probability density function for droplet radius
!
! -------------------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------------------

        SUBROUTINE radius_pdf(n_drop,pNmin,pNmax,nbin,radius,r_mean,r2_mean,r3_mean,pdf,scalar)

        IMPLICIT NONE

        INTEGER :: i, j, nbin, n_drop
        REAL*8  :: ddx, pNlength, pNmin, pNmax, pdf_area
        REAL*8  :: radius(n_drop), pdf(nbin), r(nbin), r2(nbin), r3(nbin), scalar(nbin)
        REAL*8  :: r_mean(nbin), r2_mean(nbin), r3_mean(nbin)


        DO i=1,nbin
           pdf(i) = 0.0
           r(i)   = 0.0
           r2(i)  = 0.0
           r3(i)  = 0.0
        END DO

! -- PDFAREA  the total area of PDF
! -- k        locate the A(k) into approciate bin #
! -- pNmin    minimum of scalar value
! -- pNmax    maximum of scalar value

        pNlength = pNmax-pNmin
        pdf_area  = 0.0

        DO i=1,n_drop
           j = INT(nbin*((radius(i)-pNmin)/pNlength))+1

           IF(j .GT. nbin ) THEN  
              j = nbin
           ELSE IF (j .LE. 0) THEN
              j = 1
           END IF
! -- accumulate the count # at each interval of j
           pdf(j) = pdf(j)+1.0
           r(j)   = r(j)+radius(i)*1.0e6
           r2(j)  = r2(j)+radius(i)*radius(i)*1.0e12
           r3(j)  = r3(j)+radius(i)*radius(i)*radius(i)*1.0e18
        END DO

! -- calculate the bin width
        ddx = (pNlength)/FLOAT(nbin)

! -- normlaized the PDF and calculate the real scalar value at each bin#
        DO i=1,nbin
           IF (pdf(i) .EQ. 0) THEN
              r_mean(i)  = 0.0
              r2_mean(i) = 0.0
              r3_mean(i) = 0.0
           ELSE 
              r_mean(i)  = r(i)/pdf(i)
              r2_mean(i) = r2(i)/pdf(i)
              r3_mean(i) = r3(i)/pdf(i)
           END IF
           pdf(i)    = pdf(i)/n_drop
           scalar(i) = (i-1.0)*ddx+pNmin+ddx/2.0
        END DO

        RETURN
        END

