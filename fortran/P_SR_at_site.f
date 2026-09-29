
      subroutine calc_P_SR_atSite (rupLength, fltLength, xl_F, P_SR_atSite )
      real rupLength, fltLength, xl_F, P_SR_atSite
      real P1, P2

c     Probability of the surface rupture at the site based on a uniform distribution

      if ( rupLength/fltLength . le. xl_F ) then
            P1 = (xl_F - rupLength/fltLength) / ( 1. - rupLength/fltLength)
      else
            P1 = 0.
      endif

      if ( rupLength/fltLength . le. 1.-xl_F ) then
            P2 = (1. - xl_F - rupLength/fltLength) / ( 1. - rupLength/fltLength)
      else
            P2 = 0.
      endif
      P_SR_atSite = 1. - (P1+P2)

      return
      end


