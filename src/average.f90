! -------------------------------------------------------------------------------------------------
! This subroutine calculates the averaged rms and mean of the scalars of different realizations
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

	SUBROUTINE average_2(scalar_ave,scalar_int,ip,mr,mip)

        IMPLICIT NONE

        INTEGER :: i, ip, mip, mr
	REAL*8  :: scalar_ave(8,0:mip)
	REAL*8  :: scalar_int(8,0:mip)
        
        DO i=1,8
           scalar_ave(i,ip) = scalar_ave(i,ip)*(mr-1)/mr+scalar_int(i,ip)/mr
        END DO

        RETURN

        END

! -------------------------------------------------------------------------------------------------
! This subroutine averages the stastics at realizations
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

        SUBROUTINE AVERAGE( STATave,STATint,nbin,ipe,mr,mip)

! -- STATint  the input value of instaneous value at current realization
! -- STATave  the ave. value at different real.
! -- mr       number of realization
        IMPLICIT NONE

        INTEGER :: i, ipe, mip, mr, nbin
        REAL*8  :: STATave(mip,101), STATint(mip,101)

        DO i=1,nbin!+1
           STATave(ipe,i) = STATave(ipe,i)*(mr-1)/mr+(STATint(ipe,i)/mr)
        END DO

        RETURN

        END

!c***************************************************************************
!c
!c This is a subroutine to average the scalars  at realizations
!c Note the array is increased to 2000
!c STATint: the input value of instaneous value at current realization
!c STATave: the ave. value at different real.
!c MR: number of realization
!c
!c***************************************************************************
!
!C	SUBROUTINE AVERAGE_SCALAR(STATave, STATint,M1,M4, Ipe, MR)
!C       real*8 STATave(20,2000), STATint(20,2000)
!C       DO JJ=M1,M4
!C       STATave(Ipe,JJ)=STATave(Ipe,JJ)*(MR-1)/MR+(STATint(IPe,JJ)/MR)
!C       END DO
!C       RETURN
!C       END
!
!c***************************************************************************
!C
!C       This is a subprogram to calculate Phiprime, phiprim2,
!C       phiprim4, phiprim6.
!C       Phiprime=< phi' >
!C       Phiprime2=< (phi')^2 >
!C	phiprime4=< (phi')^4 >
!C	Phiprome6=< (phi')^6 >
!C
!C***************************************************************************
!c comment out
!c
!cc       SUBROUTINE PHI(M1,M4,BB,PHIRMS,PHIPRIM2,PHIPRIM4,PHIPRIM6,igrid)
!cc	integer*4 igrid
!cc       real*8 BB(0:igrid)
!cc       G1=0.
!cc       DO I=M1,M4
!cc       G1=G1+BB(I)
!cc       END DO
!
!cc       PHIMEAN=G1/(M4-M1+1)
!cc       PHIPRIM=0.
!cc       PHIPRIM2=0.
!cc       PHIPRIM4=0.
!cc       PHIPRIM6=0.
!cc       DO I=M1,M4
!cc       PHIPRIM=PHIPRIM+(BB(I)-PHIMEAN)
!cc       PHIPRIM2=PHIPRIM2+((BB(I)-PHIMEAN))**2.
!cc       PHIPRIM4=PHIPRIM4+((BB(I)-PHIMEAN))**4.
!cc       PHIPRIM6=PHIPRIM6+((BB(i)-PHIMEAN))**6.
!cc       END DO
!
!cc       PHIPRIM=PHIPRIM/(M4-M1+1)
!cc       PHIPRIM2=PHIPRIM2/(M4-M1+1)
!cc       PHIRMS=PHIPRIM2**.5
!cc       PHIPRIM4=PHIPRIM4/(M4-M1+1)
!cc       PHIPRIM6=PHIPRIM6/(M4-M1+1)
!cc       RETURN
!cc       END
!
!****************************************************************************
!*
!*       This a subroutine to calculate the dissipation rate
!*       and the correlation.
!*
!*       BB: scalar field
!*       Depsilon: dissipation rate
!*       Dro2: correlation factor
!*       DM: diffusivity coefficent
!*
!****************************************************************************
!
!cc       SUBROUTINE DISSIPATION(M1,M4,BB,Depsilon,Dro2,DX,DM,igrid)
!ccinteger*4 igrid
!cc       real*8 BB(0:igrid)
!cc
!c 
!c Whta's DM here? Sould we change these statistics subroutine to
!c become dimensional?
!c
!c       DM=.035
!
!cc       Dro=0.
!cc       DD2=0.
!!
!
!cc       Depsilon=0.
!CC TN = FLOAT(M4-M1+1)!
!
!CC     BB(M4+1)=BB(M1)
!CC      DO 10 K=M1,M4
!c
!c       differentiate BB(K)
!c
!c       D1=(BB(K+1)-BB(K))/DX
!c       Depsilon=D1**2.+Depsilon
!c       D2=BB(K)**2.
!c      DD2=DD2+D2
!c       D3=DM*(D1**2.)
!cD3 = 1.0
!c       Dro=Dro+D2*D3
!c10      CONTINUE
!
!c       Depsilon=DM*Depsilon/TN
!c       DD2=DD2/TN
!c       Dro1=Dro/TN
!c       Dro2=Dro1/(DD2*Depsilon)-1.
!c       RETURN
!c       END
!!
!
!****************************************************************************
!*
!* This is to ave. the statistics quantities.-int is the input
!* value, -ave is the ave. value.
!* where VA1 - rms of scalar
!*       VA2 - kurtosis
!*       VA3 - superskewness
!*       VA4 - dissipation rate
!*       VA5 - correlation factor
!*       ipe _ period of interval during the realization
!*        MR _ realization
!*
!****************************************************************************
!
!
!c       SUBROUTINE AVERAGE_1(sVA1ave,VA1int,sVA2ave,VA2int
!c    1  ,sVA3ave, VA3int, sVA4ave, VA4int, sVA5ave, VA5int,ipe,MR,mip)
!
!c       DIMENSION sVA1ave(mip), VA1int(mip)
!c       DIMENSION sVA2ave(mip), VA2int(mip)
!c       DIMENSION sVA3ave(mip), VA3int(mip)
!c       DIMENSION sVA4ave(mip), VA4int(mip)
!c       DIMENSION sVA5ave(mip), VA5int(mip)
!
!
!
!c       sVA1ave(ipe)=sVA1ave(ipe)*(MR-1)/MR+VA1int(ipe)/MR
!c       sVA2ave(ipe)=sVA2ave(ipe)*(MR-1)/MR+VA2int(ipe)/(MR)
!c       sVA3ave(ipe)=sVA3ave(ipe)*(MR-1)/MR+VA3int(ipe)/(MR)
!c       sVA4ave(ipe)=sVA4ave(ipe)*(MR-1)/MR+VA4int(ipe)/(MR)
!c       sVA5ave(ipe)=sVA5ave(ipe)*(MR-1)/MR+VA5int(ipe)/(MR)
!
!c       RETURN
!c       END
!
