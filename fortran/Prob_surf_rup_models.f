

      subroutine calc_PsurRup ( Mag, P_surf_rup, SOF, jcalc_PSR, c1, c2, RWR )
      implicit none
      real mag, P_surf_rup, RWR
      integer SOF, jcalc_PSR, iModel
      real c1, c2
      character*80 name_PSR

      if (jcalc_PSR .eq. 0 ) then
        P_surf_rup =1.
      endif
      
c     User-specified logistic for Mag
      if (jcalc_PSR .eq. 1 ) then
        P_surf_rup =EXP(c1+c2*mag)/(1+EXP(c1+c2*mag))
      endif
      
c     WC PSR model for all SOF
      if (jcalc_PSR .eq. 3 ) then
        P_surf_rup =EXP(-12.51+2.053*mag)/(1+EXP(-12.51+2.053*mag))
        name_PSR = 'Wells and Coppersmith 1984'
        return
      endif
      
c     Huang and ABrahamson (2026) PSR model for all SOF
c     Median model
      if (jcalc_PSR .eq. 4 ) then
        iModel = 1
        call Huang_Abrahamson_2026_PSR ( RWR, P_surf_rup, iModel )
        name_PSR = 'Huang and Abrahamson 2026 median alpha'
        return
      endif

c     Median model
      if (jcalc_PSR .eq. 5 ) then
        iModel = 2
        call Huang_Abrahamson_2026_PSR ( RWR, P_surf_rup, iModel )
        name_PSR = 'Huang and Abrahamson 2026, 5th alpha'
        return
      endif
      
c     Median model
      if (jcalc_PSR .eq. 6 ) then
        iModel = 3
        call Huang_Abrahamson_2026_PSR ( RWR, P_surf_rup, iModel )
        name_PSR = 'Huang and Abrahamson 2026, 95th alpha'
        return
      endif

      return
      end

c ------------------------------------------

      subroutine simple_PSR ( mag, P_surf_rup )
      real mag, P_surf_rup

c     Simple test model for PSR
      if ( mag .gt. 7.5 ) then
           P_surf_rup = 1.
      elseif ( mag .ge. 5.0 ) then
           P_surf_rup = -2.88 + 0.535*mag + 0.0464*(mag-6.25)**2 -0.1*(mag-6.25)**3
      else
           P_surf_rup = 0.
      endif
      return
      end

c ------------------------------------------

      
      subroutine Huang_Abrahamson_2026_PSR ( RWR, P_surf_rup, iModel )
      real RWR, P_surf_rup
      real coeff_50(3), coeff_05(3),coeff_95(3)
      data coeff_50 /   0.2191,	1.0669,	-0.2879 /
      data coeff_05 / 0.5067,	0.7712,	-0.2802 /
      data coeff_95 /0.0590,	1.0837,	-0.1435 /
      
c     Approx to a beta distribution
      if ( iModel .eq. 1 ) then
        P_surf_rup =  coeff_50(1)*RWR + coeff_50(2)*RWR**2 +coeff_50(3)*RWR**3 
      elseif ( iMOdel .eq. 2 ) then
        P_surf_rup =  coeff_05(1)*RWR + coeff_05(2)*RWR**2 +coeff_05(3)*RWR**3 
      elseif ( iModel .eq. 3 ) then
        P_surf_rup =  coeff_95(1)*RWR + coeff_95(2)*RWR**2 +coeff_95(3)*RWR**3 
      else
        write (*,'( 2x,''Invalid PSR model index'')')
      endif

      return
      end

 
      
