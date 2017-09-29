! -------------------------------------------------------------------------------------------------
! This function calculates the saturated partial pressure of water vapor over liquid water and ice
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

        REAL*8 FUNCTION ew(tz)  

        IMPLICIT NONE

! -- ew     saturated partial pressure of water vapor over liquid water (mb)
! -- ei     saturated partial pressure of water vapor over ice (mb)
! -- tz     absolute temperature
!
! -- the formulas are from Pruppacher and Klett, p.625                     

! -- these formulas are valid from -50 C to +50 C for water and -50c to 0 c for ice.
! -- if the temperature is outside these ranges, the partial pressure is calculated at the valid
! -- temperature limit for t.gt.50C, and the partial pressure is calculated with an exponential 
! -- decay from -50 C.                  
        REAL*8 :: t, ta, tz
        REAL*8 :: a0, a1, a2, a3, a4, a5, a6, ei, et

        a0 = 6.107799961
        a1 = 4.436518521e-01
        a2 = 1.428945805e-02
        a3 = 2.650648471e-04
        a4 = 3.031240396e-06
        a5 = 2.034080948e-08
        a6 = 6.136820929e-11

        t  = tz-273.16
        ta = 1.0

        IF (t .GT. 50.0) t = 50.0

        IF (t .LT. -50.0) THEN
           ta = t
           t  = -50.0
        END IF

        ew = a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6)))))

        IF (ta .LT. -50.0) THEN
           ew = ew*dexp((ta+50D0)/10.0)
        END IF
 
        RETURN 
                                                                     
        entry ei(tz)                                                      
                                                                    
        a0 = 6.109177956                                                    
        a1 = 5.034698970e-01                                                
        a2 = 1.886013408e-02                                                
        a3 = 4.176223716e-04                                                
        a4 = 5.824720280e-06                                                
        a5 = 4.838803174e-08                                                
        a6 = 1.838826904e-10                                                
                                                                       
        t  = tz-273.16                                                       
        ta = 1.0    
                                                        
        IF (t .LT. -50.0) THEN                                               
           ta = t                                                     
           t  = -50.0                                                  
        END IF
                                                         
        IF (t .GT. 0.0) t = 0.0                                                

        et = a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6)))))     
                
        IF (ta .LT. -50.0) THEN                                              
           et = et*DEXP((ta+50.0)/10.0)                                 
        END IF                                                           
                                                                 
        ei = et    
                                                         
        RETURN
                                                            
        END                                                               
