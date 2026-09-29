      subroutine calc_fDM ( iFDM_index, Mag, iSOF, XL1, 
     1           muD, sigmaD, P_Gap, P_DP_zero, b2, mu_agg, sigma_Dagg )

      implicit none
      real Mag, XL1, mu_agg, mu_p, sigma_Dagg, sigma_DP
      real mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime
      real P_Gap, P_DP_zero
      integer iSOF, iFDM_index
      real muD, sigmaD, Dmed, sigma_NE, mu_distrib
      real b2
      
      muD = -1.E30

c     LA23 model - principal disp, single segment model
      if ( iFDM_index .eq. 1 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2 )
        muD = mu_p
        sigmaD = sigma_DP
        mu_distrib = (mu_agg**(1./0.3) - mu_p**(1/0.3))**(0.3)
        return
      endif

c     LA23 model - principal disp, unsegment model
      if ( iFDM_index .eq. 2 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2 )
        muD = mu_p_prime
        sigmaD = sigma_DP_prime
        mu_distrib = (mu_agg_prime**(1./0.3) - mu_p_prime**(1/0.3))**(0.3)
        return
      endif

c     LA23 model - aggregate disp, single segment model
      if ( iFDM_index .eq. 3 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2  )
        muD = mu_agg
        sigmaD = sigma_Dagg
        mu_distrib = (mu_agg**(1./0.3) - mu_p**(1/0.3))**(0.3)
        return
      endif

c     LA23 model - aggregate disp, unsegment model
      if ( iFDM_index .eq. 4 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2  )
        muD = mu_agg_prime
        sigmaD = sigma_Dagg_prime
        mu_distrib = (mu_agg_prime**(1./0.3) - mu_p_prime**(1/0.3))**(0.3)
        return
      endif
      
c     Non-ergodic sigma liou & abrahamson (2025) sigma_mag = 0.1
      if ( iFDM_index .eq. 5 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2  )
        Dmed = mu_p_prime**(1./0.3)
        call Liou_NE_sigma ( sigma_NE, 0.1, Dmed)
        sigmaD = sigma_NE
        return
      endif

c     Non-ergodic sigma liou & abrahamson (2025) sigma_mag = 0.15
      if ( iFDM_index .eq. 6 ) then
        call LA23_FDM ( Mag, iSOF, XL1, 
     1           mu_agg, mu_p, sigma_Dagg, sigma_DP,
     1           mu_agg_prime, mu_p_prime, sigma_Dagg_prime, sigma_DP_prime,
     3           P_Gap, P_DP_zero, b2  )
        Dmed = mu_p_prime**(1./0.3)
        call Liou_NE_sigma ( sigma_NE, 0.15, Dmed)
        sigmaD = sigma_NE
        return
      endif

        write (*,'( 2x,'' bad iFDM_index'',i5)') iFDM_index
        stop 99

      return
      end
      


    
