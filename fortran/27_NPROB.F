                                                                                                                                                      
      subroutine S27_NDTR( X, P, D)                                                 
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                                
C          D IS DENSITY VALUE (OUTPUT).

      implicit none

      real X, AX, P, D, T

      IF (X) 1,2,2                                                              
    1 AX = -X                                                                   
      GOTO 3                                                                    
    2 AX = X                                                                    
    3 IF ( AX-6.0 ) 5,4,4                                                       
    4 P = 1.                                                                    
      D = 0.                                                                    
      GOTO 6                                                                    
    5 T = 1. / (1.0 + 0.2316419 * AX)                                           
      D = 0.3989423 * EXP(-X*X / 2.0)                                           
      P = 1.0 - D*T*( (((1.330274*T - 1.821256)*T + 1.781478) * T-              
     1    0.3565638) * T + 0.3193815)                                           
    6 IF (X) 8,7,7                                                              
    7 P = 1.0 - P                                                               
    8 continue                                                                  
      RETURN                                                                    
      END                                                                                                                                             

c ------------------------------------------------------------------
      subroutine S27_NDTR3( X, P)    
C          Reference: Abramowitz and Stegan equation 7.1.26                                             
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                               

      implicit none

      real x
      real*8 p, x1, x2, p1, a1, a2, a3, a4, a5, t
      data p1, a1, a2, a3, a4, a5 / 0.3275911, 0.254829592, 
     1     -0.284496736, 1.421413741, -1.453152027, 1.061405429 /
     
      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
      endif
	  
        x2 = x1/(sqrt(2.))
        t = 1/(1+(p1*x2))
        p = 1-0.5*((a1*t)+(a2*(t**2))+(a3*(t**3))+(a4*(t**4))+(a5*(t**
     1      5)))*(exp(-(x2**2)))
     
      if ( x .gt. 0. ) then
	    p = 1. - p
      endif

      RETURN                                                                    
      END   
c ------------

      subroutine S27_NDTR3x( X, P)    
C          Reference: Abramowitz and Stegan equation 7.1.26                                             
C          X IS NO. OF STANDARDIZED NORMAL DEVIATES.                            
C          P IS COMP. CUMULATIVE VALUE (OUTPUT).                               

      implicit none

      real x
      real*8 p, x1, x2, p1, a1, a2, a3, a4, a5, t
      data p1, a1, a2, a3, a4, a5 / 0.3275911, 0.254829592, 
     1     -0.284496736, 1.421413741, -1.453152027, 1.061405429 /
     
      if ( x .lt. 0. ) then
	    x1 = abs(x)
	  else
	    x1 = x
      endif
	  
        x2 = x1/(sqrt(2.))
        t = 1/(1+(p1*x2))
        p = 1-0.5*((a1*t)+(a2*(t**2))+(a3*(t**3))+(a4*(t**4))+(a5*(t**
     1      5)))*(exp(-(x2**2)))
c        write (*,'( e12.3)') a1, a2, a3, a4, a5, t, x2, p
c        pause 'in nprob'
     
      if ( x .gt. 0. ) then
	    p = 1. - p
      endif

      RETURN                                                                    
      END       
    
c ------------------------------------------------------------------

      real function pxceed ( eti, ti, siga, sigTrunc )                   

      implicit none
            
      integer jj, i                                
      real eti, ti, siga, sigTrunc, W, G, D, g1 
                                                                                      
      W = (ti-ETI)/SIGA                                                   
      call S27_NDTR(W,G,D)                                                          
                                                                                
c     TRUNCATE DISTRIBUTION AT SigTrunc                                         
      if (w .gt. sigTrunc) then                                                 
        pxceed = 0.0                                                                   
      else                                                                      
        call S27_NDTR(sigTrunc,g1,D)                                                       
        pxceed = (g-g1)*(1.+g1)                                                        
      endif 
      return                                                                    
      end                                                                        

c ------------------------------------------------------------------

      real*8 function pxceed3 ( eti, ti, siga, sigTrunc ) 

      implicit none                    
       
      integer jj, i                                           
      real eti, ti, siga, sigTrunc, w
      real*8 g, g1

c     truncates and renormalizes on both low and high end
                                                                                
      w = (ti-eti)/siga                                                   
      call S27_NDTR3(w,g)     
                                                                                                                         
      if (w .gt. sigTrunc) then                                                 
        pxceed3 = 0.0  
      else if (w .lt. (sigTrunc*(-1.))) then
        pxceed3 = 1.0                                                                 
      else                                                                      
        call S27_NDTR3(sigTrunc,g1)                                                       
        pxceed3 = (g-g1)/(1.-(2.*g1))                                                             
      endif 
                                                                    
      return                                                                    
      end                                                                                                                                             

c ------------------------------------------------------------------

 
