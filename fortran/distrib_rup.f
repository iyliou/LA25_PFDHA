 
c      ------------------------------------------------------------------

      subroutine calc_pXceed_DD ( iSOF1, DTot_rat, pxceed_DD ) 
      real DTot_rat, pxceed_DD
      integer iSOF, ISOF1
      real d4(3), d5(3), d6(3)
      data d4 / -0.855, -0.968, -0.362 /
      data d5 / 0.262, 0.403, -0.159 /
      data d6 / 0.076, 0.018, 0.194 /
      
      iSOF = iSOF1 + 2
      Pxceed_DD = 1 + d4(ISOF)*DTot_Rat + d5(ISOF)*DTot_Rat**2 
     1     + d6(ISOF)*DTot_Rat**3 
      return
      end

c -----------------------------------------------

      subroutine calc_P_atsite_W ( Rx1, F_map, W_site, iSOF, P )
      implicit none
      real Rx, w_site, P, Rx1
      real eps1, eps2, fx_RX
      real d2(3), d3(3)
      integer ISOF, F_map, ISOF1
      real*8 P1, P2
      data d2 / 0.976, 0.742, 0.942 /
      data d3 / 761, 1515, 647 /

c     Rx1 is in km
c     W_site is in m

      ISOF1 = ISOF + 2

      Rx = Rx1 * 1000.

      eps1 = (rx + W_site) / d3(ISOF1)
      eps2 = Rx / d3(ISOF1)
      call S27_NDTR3x( eps1, P1)
      call S27_NDTR3x( eps2, P2)
c      write (*,'( 2e10.3)') eps1, P1
c      write (*,'( 2e10.3)') eps2, P2
      if ( Rx1 .lt. 15. ) then
       fx_RX = 1/17.5
      elseif ( Rx1 .lt. 20. ) then
       fx_RX = (1. - (Rx1-15.)/ 5.)/17.5
      else
        fx_RX = 0.
      endif
      
      P = 2.*d2(ISOF1) *(P2-P1) + (1.-d2(ISOF1))*F_map * fx_RX * W_site/1000.
c      write (*,'(e10.4)') d2(ISOF1), P1, P2, (1.-d2(ISOF1))
c      write (*,'(e10.4)') fx_RX , W_site, P
c      pause 

      return
      end

c -------------------------------------------------------

      subroutine calc_D_DTot (P_D_tot, D_Dtot, eps_step, b2_eps, iSOF, Dagg )
      implicit none
      real P_D_tot, D_DTot, eps_step, b2_eps, Dagg
      real D_Dtot_PN
      integer ISOF, iSOF1
      real b2(3), phi_b2, tau_b2, sigma_b2
      real*8 P1, P2, psum
      real muP_prime, mu_D_Dtot, eps1, eps2
      
      data b2 / -0.097, -0.048, -0.062/
      
      common /temp1/ psum
      
      ISOF1 = ISOF + 2
      
      phi_b2 = 0.11
      tau_b2 = 0.05
      sigma_b2 = sqrt(phi_b2**2 + tau_b2**2)
            
c     set the median total distributed disp
      muP_prime =  max ( Dagg**(0.3) + b2(ISOF1), 0.)
      mu_D_DTot = Dagg**(0.3) - muP_Prime
c      write (*,'( ''Dagg, muP_prime, mu_D_Dtot'', 3f10.3)') Dagg, muP_prime, mu_D_Dtot
      
c     compute D_Dtot_PN (in m^0.3 units) for this epsilon
      D_Dtot_PN = ( mu_D_DTot + sigma_b2 * b2_eps) 
c      write (*,'( 3f10.4)') D_Dtot_PN, b2_eps, sigma_b2

c     Convert back to linear units
      if ( D_Dtot_PN .le. 0. ) then
        D_Dtot =0.
      else 
        D_Dtot = D_Dtot_PN**(1./0.3)      
      endif
c      write (*,'( e12.3)') D_Dtot

c     D_Dtot can't be greater than D_agg. This truncates the distribution
      if ( D_Dtot .gt. Dagg ) D_Dtot = Dagg
c      write (*,'( ''after truncation'',e12.3)') D_Dtot

c     compute probability of mu_Dtot between (x-step/2 and x+step/2)
      eps1 = b2_eps - eps_step/2.
      eps2= b2_eps + eps_step/2.
      call S27_NDTR3x( eps1, P1)
      call S27_NDTR3x( eps2, P2)
      P_D_tot = P1-P2
c      write (*,'( 2x,''P1, P2, P'',10f10.4)') b2_eps, eps_step, eps1, eps2, P1, P2, P_D_tot
c      pause
          
      Psum = psum + P_D_tot
      
      return
      end

c --------------------

      subroutine PD_petersen ( Rx, L_site, W_site,   PD_GT_0_atSite )   
      implicit none
      real Rx, L_site, W_site,   PD_GT_0_atSite 
      real P0, P1, P2, PD, a, b, sum1, Rx1
      integer i, N1
      
      P0 = 0.902
      P1 = 0.185
      P2 = 0.0664

c     Petersen model for 100m x 100 m
      a = -1.0114
      b = 2.557
      
      sum1 = 0.
      N1 = int( W_site / 10. )
      if ( N1 .eq. 0 ) N1 = 1
      do i=1,N1
        Rx1 = Rx + (i-1)*10.
        if ( Rx .lt. 0.1 ) then
          PD = Rx/0.1 * (P1-P0) + P0
        elseif ( Rx .lt. 0.2 ) then
          PD = (Rx-0.1)/0.1 * (P2-P1) + P1
        else
          PD = exp(a*log(Rx*1000.) + b)
        endif
        sum1 = sum1 + PD/10.
      enddo
      
c     adjust prob for the length of the site, assuming constant along strike      
      PD_GT_0_atSite = sum1 * L_site / 100.

      return
      end
      

        
