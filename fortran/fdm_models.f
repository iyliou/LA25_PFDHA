      subroutine LA23_fDM ( Mag, iSOF1, X_over_L, 
     1       mu_agg, mu_p, sigma_dagg, sigma_Dp,
     1       mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3       P_Gap, P_DP_zero, b2_1 )

      implicit none
      real c0, c1, c2, c3, c4, c7
      real c5(3), c6(3), c10(3), c11(3), c12(3), c13(3), c14(3), 
     1     c14a(3), c15(3), c16(3), c17(3), c17a(3),
     1     c18(3), c19(3), c20(3), c21(3), c22(3), c23(3),
     2     c24a(3), c24b(3), c24c(3), c25(3), c26a(3), c26b(3), c26c(3)
      real M1(3)
      real b0(3), b1(3), b2(3), b3(3), b4(3), b5(3), phi_b2
      real e1(3), e2(3), e3(3)
      real*8 pxceed3
      real Mag, X_over_L, XL_1, XL, z, Pz
      real phi, phi1, phi2, tau1, tau2, tau, sigTrunc
      real mu_agg, D_K_agg, T_M, T_XL
      real c24_mag, c26_mag
      real f_N_Delta_Mu, deltaMu_max,  mu_agg_prime, sigma_add
      real P_GAP_max,  f_N_P_gap, P_GAP, P_DP_zero, mu_P_prime, mu_p
      integer iSOF, iSOF1, iRefSOF(3), k
      real Dagg_med, Dp_med, sigma_Dagg, sigma_DP, b2_med, b2_sigma
      real  Dagg_med_prime, Dp_med_prime, sigma_DAgg_prime, sigma_Dp_prime
      real b2_1

      data iRefSOF / 1, 2, 3 /
      data c5 / 0.7, 0.7, -0.2 /
      data c6 / -0.262, -0.314, 0.0 /
      data M1 / 6.25, 6.5, 5.0 /
      data c10 / -0.406, -0.394, -0.849 /
      data c11 / -8.205, -2.95, 8.026 /
      data c12 / -165, -64.3, 237.2 /
      data c13 / -1372, -833, 927 /
      data c14 / -3013, -2145, 831 /
      data c14a / -0.904, -1.01, -3.555 /
      data c15 / -0.44, -0.155, -0.087 /
      data c16 / 0.0767, 0.0263, 0.0148  /
      data c17 / 0.0275, 0.0177, 0.0087 /
      data c17a / 0.00276, 0.00534, 0.00234 /
      data c18 / -0.215, -0.201, -0.165 /
      data c19 / 0.043, 0.0364, 0.0304 /
      data c20 / 0.00574, 0.0102, 0.00773 /
      data c21 / -0.082, -0.135, -0.151 /
      data c22 / 0.027, 0.027, 0.033 /
      data c23 / -0.0088, 0.0063, 0.0032/
      data c24a / 0.02, 0.04, 0.04 /
      data c24b / 0.00, 0.005, 0.00 /
      data c24c / 0.0247, 0.00273, 0.00843 /
      data c25 / 0.6, 0.31, 0.43 /
      data c26a / 0.162, 0.153, 0.156 /
      data c26b / 0.084, 0.017, 0.031 /
      data c26c / 0.008, 0.00312, 0.00575 /

      data  b0 /-3.24,0.867, 1.652/
      data  b1 /9.105, 3.767, 1.349/
      data  b2 /-0.097, -0.048, -0.062/
      data  b3 /0.286, 0.191, 0.28 /
      data  b4 /1.195, 1.058, 0.0 /
      data  b5 /-1.526, -1.418, 0.0/
      
c      data  e1 /0.815, 1.216, 0.646/
c      data  e2 /0.354, 0.177, 0.432/
c      data  e3 /-0.0958, -0.0511, -0.115/    

      iSOF = iRefSOF(iSOF1+2)

c     Set coeff that are independent of SOF 
      c0 = 0.272
      c1 = 0.913
      c2 = -2.128
      c3 = 1.15
      c4 = 0.0
      c7 = 0.033
      tau1 = 0.115
      tau2 = 0.205
      phi1 = 0.12
      phi2 = 0.27
      phi_b2 = 0.11
      
c     fold X_over_L to [0,0,5] range
      XL = X_over_L
      if ( XL .gt. 0.5 ) then
        XL = 1-X_over_L
      endif

c     Compute the med displacement from the wavenumeber model before tapers (eq 8 in LA23)
      D_K_agg = exp(c0 + c1*(XL-0.3) + c2*(XL-0.3)**2 + c3*(Mag-7.))
c      write (*,'( ''D_K_agg='',f10.4)') D_K_agg


c     Compute XL1 (eq 10 in LA23)
      if ( mag .lt. 7. ) then
        XL_1 = 0.15
      elseif ( mag .lt. 8. ) then
        XL_1 = 0.15 - 0.1*(mag-7.)
      else
        XL_1 = 0.05
      endif
        
c     Compute X/L taper (eq 9 in LA23)
      if ( XL .lt. XL_1 ) then
        T_XL = (XL - XL_1) /XL_1
      else
        T_XL = 0.
      endif


c     Compute M taper (eq 13 in LA23)
      if ( mag .gt. 7.) then
        T_M = 0.
      elseif ( mag .gt. M1(iSOF) ) then
        T_M = (7.-mag)/(7.-M1(iSOF))
      else
        T_M = 1.0
      endif

c     compute med Disp from  tapered K model on m^0.3 units (eq 14)
      mu_agg = ( D_K_agg * exp(c5(iSOF)*T_XL))**(0.3) + c6(iSOF)*T_M + c7
c      write (*,'( '' mu_agg='',3f10.4)') XL, T_XL,  mu_agg
c      pause
      
c     compute Phi (eq 15)
      if ( mag .lt. 6 ) then
        phi = 0.12
      elseif (mag .lt. 7 ) then
        phi = 0.12 + 0.15*(Mag-6.)
      else
        phi = 0.27
      endif

c     compute Tau (eq 16)
      if ( mag .lt. 6 ) then
        tau = 0.115
      elseif (mag .lt. 7.5 ) then
        tau = 0.115 + 0.06*(Mag-6.)
      else
        tau = 0.205
      endif
            
c     compute adjustment to median (delta_mu_agg) to account for segments
c     compute normalized X/L shape for deltaMu (eq 20)
      if ( XL .le. 0.3 ) then
        f_N_Delta_Mu = c10(iSOF) + c11(iSOF)*(XL-0.3) + c12(iSOF)*(XL-0.3)**2
     1         + c13(iSOF) *(XL-0.3)**3 + c14(iSOF)*(XL-0.3)**4
      elseif ( XL .le. 0.4) then
        f_N_Delta_Mu = c10(iSOF)
      else
        f_N_Delta_Mu = c10(iSOF) + c14a(iSOF)*(XL-0.4)
      endif
      
c     compute deltaMu_max (eq 21)
      deltaMu_max = c15(iSOF) +c16(iSOF)*mag + c17(iSOF)*(mag-6.7)**2
     1             + c17a(iSOF)*(Mag-6.7)**3 
     
c     COmpute the mu_agg' (equation 22)
      mu_agg_prime = mu_agg + deltaMu_max * f_N_Delta_Mu  
c      write (*,'( 2f10.4)') deltaMu_max , f_N_Delta_Mu  
c      write (*,'( 4f10.4)') c15(iSOF) ,c16(iSOF)*mag , c17(iSOF)*(mag-6.7)**2
c     1             , c17a(iSOF)*(Mag-6.7)**3 

c      write (*,'( ''deltaMu_max, f_N_Delta_Mu'',3f10.4)') deltaMu_max,
c     1       f_N_Delta_Mu, deltaMu_max * f_N_Delta_Mu 
c      write (*,'( '' mu_agg_prime='',f10.4)')  mu_agg_prime

      sigma_dagg = sqrt( phi**2 + tau**2)
      sigma_dp = sqrt( sigma_dagg**2 + phi_b2**2 - 0.3*phi_b2*sigma_dagg)
      
c     Compute the adjusted sigma (eq 17)
c     compute sigma_add (equation 17)
      sigma_add = c18(iSOF) + c19(iSOF)*mag + c20(iSOF)*(Mag-6.7)**2
c      write (*,'(i5,3f10.4)') isof, c18(isof), c19(isof), c20(isof)
      
c      write (*,'( 4f10.3)') c18(iSOF) + c19(iSOF)*mag, 
c     1      c20(iSOF)*(Mag-6.7)**2
c      write (*,'( 4f10.4)') phi, tau, sigma_dagg, sigma_add

c     compute sigma_agg prime (eq 18)
      sigma_Dagg_prime = sqrt( phi**2 + tau**2 + sigma_add**2)
c      write (*,'( 2x,''sigma agg prime'',f10.4)') sigma_agg_prime
      sigma_dP_prime = sqrt( sigma_Dagg_prime**2 + phi_b2**2 -0.3*phi_b2*sigma_Dagg_prime)

c     Compute probability of being in a gap between segments
c     compute P(GAP)_max (eq 26)
      P_GAP_max = c21(iSOF) + c22(iSOF)*mag + c23(iSOF)*(mag-6.5)**2      

c     Compute 24 and c26 (eq 28)
      c24_mag = c24a(iSOF) + c24b(iSOF)*(mag-5.) + c24c(iSOF)*(mag-5)**3
      c26_mag = c26a(iSOF) + c26b(iSOF)*(mag-5.) + c26c(iSOF)*(mag-5)**3
      
c     compute normalized X/L shape: fN_PGap (eq 27) 
      if ( XL .lt. 0.1 ) then
        f_N_P_gap = 0.
      elseif ( XL .lt. 0.15) then
        f_N_P_gap = 20. * c24_mag*(XL-0.1)
      elseif ( XL .lt. 0.28) then
        f_N_P_gap = c24_mag + 7.69*(1.-c24_mag)*(XL-0.15)
      elseif ( XL  .lt. 0.3) then
        f_N_P_gap = 1.0
      elseif ( XL .lt. 0.4) then
        f_N_P_gap = 1. - 10. *(1-c25(iSOF))*(XL-0.3)
      else
        f_N_P_gap = c25(iSOF) - 10.*(c25(iSOF)-c26_mag) *(XL-0.4)
      endif
    
c     compute P_GAP (eq 25)      
      P_GAP = P_GAP_max * f_N_P_gap 
c      write (*,'( 2x,''P_Gap'',f10.4)') P_GAP

c     Compute probability of zero displacement along a segment, but not in
c     a gap between segments (eq 32)
      P_DP_zero = 1. / ( 1.+exp(b0(iSOF)+b1(iSOF)*mu_agg))
c      write (*,'( 2x,''P_DP_zero'',f10.4)') P_DP_zero

c     compute the median principal displacement if non-zero (eq. 33)
      if ( mu_agg .gt. b2(iSOF)) then
        mu_P = mu_agg + b2(iSOF)
        mu_P_prime = mu_agg_prime + b2(iSOF)
      else
        mu_P_prime = 0.
      endif
c      write (*,'( 2x,''mu_P_prime'',f10.4)') mu_P_prime

      b2_1 = b2(iSOF)
      
      return
      end
      
c ----------------------------------------------------------------

      subroutine Liou_NE_sigma ( sigma_NE, mag_sigma, Dmed)
      implicit none
      real sigma_NE, mag_sigma, Dmed
      real sigma_NE_MX, Dmean, ratio, sigma_MX
      real d1, c0, c2


c     compute sigma_NE_MX (feb 25 model)
c       - convert median disp to mean disp (m) using the sigma_NE,MX
      ratio = 1.12  -0.10 * alog(Dmed) + 0.0223*(alog(Dmed)+0.375)**2
      Dmean = dMed * ratio
      if ( Dmean .lt.  1.44 ) then
        sigma_NE_MX = 0.17
      else
        sigma_NE_MX = 0.17 + 0.0325*(alog(Dmean/1.44))**2
      endif

      if ( dMed .lt. 0.05 ) dMed = 0.05
      
c     compute sigma_MX
      if ( mag_sigma .eq. 0.1) then
        d1 = 0.38
        c0 = 0.087
        c2 = -0.0178
        if ( dMed .gt. d1 ) then
          sigma_MX = c0
        else
          sigma_MX = c0 + c2*(alog(dMed/d1))**2
        endif
      elseif ( mag_sigma .eq. 0.15 ) then
        d1 = 0.383
        c0 = 0.0994
        c2 = -0.0206
        if ( dMed .gt. d1 ) then
          sigma_MX = c0
        else
          sigma_MX = c0 + c2*(alog(dMed/d1))**2
        endif
      else
        write (*,'( 2x,'' invalid mag_sigma'')')
        stop 99
      endif

c      write (*,'( f10.3)')  sigma_NE_MX, dMed,  sigma_MX
      
      
      if ( sigma_MX .gt. sigma_NE_MX ) then
        write (*,'( 3f10.4)') sigma_MX, sigma_NE_MX, dMed
        stop 99
      endif 
      

c     Compute sigma_NE
      sigma_NE = sqrt( sigma_NE_MX**2 - sigma_MX**2 )
      return
      end

      
      
 
