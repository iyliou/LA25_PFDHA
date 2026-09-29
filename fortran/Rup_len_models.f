c -------------------------------------------------------

      subroutine calc_RupLength ( Mag, thick, dip, SOF, P_surf_rup, jcalc_RLM,
     1                rupLength_med, sigma_RL, nRL, eps_RL0, fixedLength )
      implicit none
      real mag, dip, thick
      real F_DS, log10_Len, rupLength_med, sigma_RL, wLimit
      real P_surf_rup,fixedLength, eps_RL0
      integer SOF, jcalc_RLM, nRL

      if ( jcalc_RLM .ne. 0 ) then
        eps_RL0 = -3.
        nRL = 13
      else
        rupLength_med = fixedLength
        eps_RL0 = 0.
        nRL= 1
        sigma_RL = 0.001
        return
      endif

c     Lavrentidias et al 2022
      if (jcalc_RLM .eq. 1 ) then

c       Set down-dip width limit
        wLimit = thick / sin(dip*3.1415926/180.)

c       Set Dip-slip flag (1 for DS and 0 for SS)
        if ( SOF .eq. 0 ) then
         F_DS = 0.
        else
         F_DS = 1.
        endif 
        log10_Len = -3.69+0.92*(mag-MIN(-1.16+0.34*mag,alog10(wLimit)))-0.157*F_DS
        RupLength_med = 10**(log10_Len)
        sigma_RL = 0.26
        return
      endif

      write (*,'( 2x,'' invalid index for Rupture Length Model'',i5)')  jcalc_RLM
      stop 99
      end

c --------------
