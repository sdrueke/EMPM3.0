!C***********************************************************************
!C
!C       FUNCTION RAND1(IX)
!C
!C       This is the heart of all parts of the program that have anything
!C       to do with random numbers. It was developed by L. Schrage and
!C       taken from "Simulation Modeling and Analysis" by Averill M. Law.
!C
!C       VARIABLES:
!C
!C       IX  : Originally set to the initial seed and returns the next
!C           integer, ZSUBI, which is the seed for the next call to
!C           function RAND.
!C       RAND: Contains the next value of USUBI, which is the uniform
!C           random number being requested on the interval [0,1].
!C
!C***********************************************************************
!C
        REAL*8 FUNCTION RAND1(IX)
        IMPLICIT NONE
        INTEGER*4 A,B15,B16,FHI,IX,K,LEFTLO,P,XALO,XHI
        DATA A/16807/, B15/32768/, B16/65536/, P/2147483647/
        XHI = IX / B16
        XALO = (IX - XHI * B16) * A
        LEFTLO = XALO / B16
        FHI = XHI * A + LEFTLO
        K = FHI / B15
        IX = (((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16) + K
        IF(IX.LT.0) IX = IX + P
        RAND1 = FLOAT(IX) * 4.656612875E-10
        RETURN
        END
