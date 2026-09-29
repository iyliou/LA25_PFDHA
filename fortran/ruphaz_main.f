      program Fault_rup_hazard

      implicit none
      include "rupHaz.inc"

c     dimensions for inputs
      real rate(MAX_NMAG,MAX_BR,MAXWIDTH), rate_wt(MAX_BR,MAXWIDTH)
      real mag1(MAXWIDTH,MAX_NMAG)
      real fltLength(MAXWIDTH), dip(MAXWIDTH), seismoThick(MAXWIDTH), flt_width_DD(MAXWIDTH)
      real PSR_wt(MAX_PSR), RLM_wt(MAX_RLM), SOF_wt(3), P_loc
      real z(MAXZ), XL_F, width_wt(MAXWIDTH)
      real MagbinStep, EpsBinStep, mag0, eps0, XLBinStep, XL0
      real deltaMU(MAX_FDM), delta_sig(MAX_FDM), FDM_wt(MAX_FDM)
      real fixedLength, rupLength
      real c1(MAX_PSR), c2(MAX_PSR), t1, t2
      real PSR_save(MAX_PSR,MAX_NMAG)   
      real muDX, P_GapX, P_DP_zeroX
      real fix_sigma(MAX_FDM)
      
      real Rup_area, AspectRatio,  Rup_width, RWR 


      integer jcalc_FDM(MAX_FDM), jcalc_Sig(MAX_FDM)
      integer iMag, nMag(100), nMag2
      integer iSOF1(3), iSOF
      integer nSOF, nBr(MAXWIDTH), iBr
      integer nWidth, iWidth
      integer nMagBin, nEpsBin, nXLBin
      integer iPSR, iRLM, nPSR, nRLM, IRL, nRL
      integer jcalc_RLM(MAX_RLM), jcalc_PSR(MAX_PSR)
      integer iFDM, iMagBin, iEpsBin, ixlBin, jSOF
      integer i, iZ, nZ, iLoc
      integer iBR_total, nEndBr

c     distributed rupture param
      real zd(maxZ), Rx, W_site, L_site
      integer F_map, iFlag_DD_Prob
      integer  nFDM, nZd, iZd
      real FP, DT_Prim, DT_Distrib
      real P_distributed, P_Rd
      real  site_width, PD_site1, PD_site2, DAMP_ratio, P_Ratio_Damp
      integer nRd_max, nRd, nRatio_distrib, iRatioD, i1
      real D2_med, DT_Distrib1, sum2
      real phi, tau
      real alpha, beta, sigma_DT, CR_DTmed_Prim
      real  DT_over_DTMean, Prob_FP
      real b2, mu_Dagg

c     hazard declaraions      
      real*8 Haz(MAXZ), marginalHaz, Mhaz(MAX_NMAG,MAXZ), Haz1(MAXZ,Max_TotalBr)
      real*8 marginalHaz_a, margHaz1
      real*8 Haz_mag(MAXZ,MAX_NMAG)
      real*8 haz_noGap_prob(MAXZ)
      real*8 Haz1D(MAXZ,Max_TotalBr)
      real*8 sum_b2(MAX_NMAG), sum_b3(MAX_NMAG), sum_b4(MAX_NMAG)
      real*8 marginalHaz_c, marginalHaz_d
      real*8 Pxceed3, rate0, sum0, rate1, rate2, rate3, rate4
      real endBR_wt(Max_TotalBr), endBR_wt4(Max_TotalBr)
      real sumwt
      real*8 Haz_FDM(MAXZ,MAX_FDM), z_med

      real*8 HazS(2,MAXZ), HazS1(2,MAXZ)
      real*8 sum33


      real meanMag(MAXZ), meanEps(MAXZ), meanXL(MAXZ)
      real X_over_l, XL1, x_site0
      real sigTrunc
      real Prob_RL
      real RupLength_med, eps_RL, sigma_RL, P13(13)
      real sigma_Dagg, P_Gap, P_DP_zero
      real mu_agg, mu_p, sigma_Dp, mu_agg_prime, mu_p_prime
      real sigma_Dagg_prime, sigma_DP_prime
      real ti,  Pz, epsilonD
      real*8 meanhaz(MAXZ)

      real P_surf_rup ,mag
      real muD, sigmaD
      
      real dmag1, minmag1, maxmag1, marg_disag_Mag(MAXZ,100)
      real deps1, mineps1, marg_disag_eps(MAXZ,100)
      real dXL1, minXL1, marg_disag_XL(MAXZ,100)
      integer nMag1, nEps1, nXL1, iMagBin1, iEpsBin1, iXLBin1
      integer nLoc2a, nLoc2b
      real x_site0a, x_site0b

c     distrib rup param
      real*8 psum, hazD(maxZ)
      real dagg_center(maxZ), dagg_rate, D_Dtot, PD_GT_0_atSite, P_D_tot
      real Dagg_rate0
      real  Pxceed_DD, b2_eps_step, b2_eps
      real deagg1(MAXZ,MaxMagBin,MaxEpsBin), eps
      real deagg2(MAXZ,MaxMagBin,MaxXLBin)
      real haz_frac(MAXZ,20)
      integer nLoc, nLoc2, ii, nAgg, iAgg
      real P_SR_atSite, x1, x2, dx
      real H1(4), HP(4)
      real LT, Prob_LT, LT_ratio, P_LT_Ratio
      real P_Li_ratio, P_order, P_SegLoc
      real P_rupLoc, length1, P1, bias, P2
      real D_Distrib_total, D2
      real epsD, eps0D, eps1D, eps2D, deps
      real D_agg1, Prob_D_agg

      data P13 / 0.0024, 0.0092, 0.0278, 0.0656, 0.121, 0.175, 0.198,
     1            0.175, 0.121, 0.0656, 0.0278, 0.0092, 0.0024 /
                                
                
c     Read inputs
      call rd_input ( nFDM, jcalc_FDM, jcalc_Sig, deltaMU, delta_sig, FDM_wt,
     1       nPSR, jcalc_PSR, PSR_wt, nRLM, jcalc_RLM, RLM_wt,
     2       nWidth, width_wt, nMag, mag1, nBR, rate_wt, rate,
     3       nSOF, iSOF1, SOF_wt, fltLength, seismoThick, dip, flt_width_DD,
     4       nZ, z, xl_F, Rx, w_site, L_site, F_map, iFlag_DD_prob, nZD, zD, 
     5       magBinStep,nMagBin, mag0, epsBinStep, nEpsBin, eps0,
     6       XLBinStep, nXLBin, XL0, fix_sigma, c1, c2, maxmag1, minMag1)
      write (*,'( 2x,''end of read input'')')

c     Initialize hazard arrays
c     Haz() is mean hazard for Dp if Rx=0 and for Dagg if Rx != 0
c     HazD() is the mean hazard for Distributed displacement
c     The meanMag, meanEps, and meanXL are for the Dagg or the Dp depending on Rx
c     Haz1() is the Dagg or Dp hazard by epistemic branch (for fractiles), depending on Rx
c     Haz1D() is the Distirbuted Disp hazard by epistemic branch (for fractile
      do iz=1,nZ
        Haz(iz) = 0.
        HazD(iz) = 0.
        HazS(1,iz) = 0.
        HazS(2,iz) = 0.
        meanMag(iz) = 0.
        meanEps(iz) = 0.
        meanXL(iz) = 0.
        
c       Initialize marginal hazard for each mag
        do iMag=1,MAX_NMAG
          MHaz(iMag,iz) = 0.
        enddo

c       initialize disagg joint
        do imagBin=1,nMagBin
         do iepsBin=1,nEpsBin
           deagg1(iz,iMagBin,iEpsBin) = 0.
         enddo
         do iXLBin=1,nXLBin
           deagg2(iz,iMagBin,iXLBin) = 0.
         enddo
        enddo

c       initiaize marginal disagg
        nMag1 = 0
        do iWidth=1,nWidth
          nMag1 = max(nMag1,nMag(iWidth))
        enddo

        deps1 = 0.1
        mineps1 = -3.
        neps1 = 61
        do i=1,neps1
          marg_Disag_Eps(iz,i) = 0.
        enddo

        dXL1 = 0.01
        minXL1 = 0.
        nXL1 = 51
        do i=1,nXL1
          marg_Disag_XL(iz,i) = 0.
        enddo

c       init second haz
        do i=1,Max_TotalBr
          Haz1(iz,i) = 0.
          Haz1D(iz,i) = 0.
        enddo
        
c       initialize haz by FMD
        do iFDM=1,NFDM
          Haz_FDM(iz,iFDM) = 0.
        enddo
        
        
      enddo
      
c     Initialize sum for hazard w/o gap prob
      do iz=1,nZ
        haz_noGap_prob(iZ) = 0.
      enddo

c     Iniitalize sums for prob of surface rup at site by Mag
      do iMag=1,MAX_NMAG
         sum_b2(iMag) = 0.
         sum_b3(iMag) = 0.
         sum_b4(iMag) = 0.
      enddo

c     write header
      write (30,'( 2x,''---------------------------------------------'')')
      write (60,'( 2x,''---------------------------------------------'')')
      if ( Rx .eq. 0. ) then
        write (30,'( 2x,''Principal displacement hazard (m) , Rx=0'')') 
        write (60,'( 2x,''Principal displacement hazard (m) , Rx=0'')') 
      else
        write (30,'( 2x,''Aggregate displacement hazard(m), Rx='',f8.3)') Rx
        write (60,'( 2x,''Aggregate displacement hazard(m), Rx='',f8.3)') Rx
      endif
      write (30,'( 2x,''---------------------------------------------'',/)')
      write (60,'( 2x,''---------------------------------------------'',/)')

      write (60,'( 18x,''mag'', 100f12.4)') (z(iz),iz=1,nz)

      iBR_total = 0
      sum0 = 0.
            
c     Loop over alternative FDMs (epistemic)
      do iFDM=1,nFDM
        write (*,'( 2x,''check sig'',2i5,f10.3)') iFDM, jcalc_sig(iFDM), fix_sigma(iFDM)

       write (*,'( 2x,''computing hazard for FMD'',i5,'' of '', i5)') iFDM, nFDM 

c      Loop over alternative PSR models (epistemic)
       do IPSR=1,nPSR

c       Loop over alternative Rup Length models (epistemic)
        do iRLM=1,nRLM

c        Loop over SOF (epistemic)
         do jSOF=1,nSOF
          iSOF = iSOF1(jSOF)

c         loop over seismogenic width (epistemic)             
          do iWidth=1,nWidth
c          set the site in km from the end of the fault             
           x_site0 = fltLength(iWidth) * xl_F
           x_site0a = fltLength(iWidth) *0.1
           x_site0b = fltLength(iWidth) * 0.9

            
c          loop over logic tree branches for the rate (epistemic)
           do iBR=1,nBR(iWidth)
           
c           set index over all logic trees and set end branch wt
c           Rate4 does not include the FDM wt. It is for computing the med disp given surface rupture by FDM
            iBR_total = iBR_total + 1                
            endBR_wt(iBR_total) = FDM_wt(iFDM) * PSR_wt(iPSR) * RLM_wt(IRLM)
     1            * SOF_wt(jSOF) * width_wt(iWidth) * rate_wt(iBR,iWidth)
            endBR_wt4(iBR_total) = PSR_wt(iPSR) * RLM_wt(IRLM)
     1            * SOF_wt(jSOF) * width_wt(iWidth) * rate_wt(iBR,iWidth)
       
c           Loop over mag (aleatory)
            do iMag=1,nMag(iWidth)
             mag = mag1(iWidth,iMag)

             sum0 = sum0 + rate(iMag,iBR,iWidth) * endBR_wt(iBR_total)

c            set the mag bin for disagregation
             iMagBin = int((mag-mag0)/MagBinStep ) + 1 
             if (iMagBin .lt. 1 ) iMagBin=1
             if (iMagBin .gt. nMagBin ) iMagBin=nMagBin

c            Compute probability of surface rupture
c            Temporary use simple model for rupture width for computing the RWR
             Rup_area = 10.**(Mag-4) 
             AspectRatio = 1.5
             Rup_width = Rup_area / AspectRatio 
             RWR = Rup_width / flt_width_DD(iWidth)
             if ( RWR .gt. 1. ) RWR = 1
             Call calc_PsurRup ( Mag, P_surf_rup, iSOF, jcalc_PSR, c1(iPSR), c2(iPSR), RWR)
             PSR_save(iPSR,iMag) = P_surf_rup
        
c            compute median and sigma of rup length
             call calc_RupLength ( Mag, seismoThick(iWidth), dip(iWidth), iSOF, P_surf_rup, jcalc_RLM(iRLM),
     1                rupLength_med, sigma_RL, nRL, eps_RL, fixedLength )

c            Integrate over rupture length (aleatory)
             do iRL=1,nRL
        
c             set rupture length, but limited to the fault length
              rupLength = 10**(alog10(rupLength_med) + eps_RL*sigma_RL) 
              if ( rupLength .gt. fltLength(iWidth)) rupLength=fltLength(iWidth)

c             Set probability of rupture length (this is for fixed sampling of eps_RL)
              if ( jcalc_RLM(iRLM) .eq. 0 ) then
               Prob_RL = 1.
              else
               Prob_RL = P13(iRL)
              endif

c             Increment epsilon for next length
              eps_RL =  eps_RL + 0.5
         
c             Compute probability of the surface rupture going past the site
c             assumes uniform distribution of the location along strike
              call calc_P_SR_atSite ( rupLength, fltLength(iWidth), xl_F, P_SR_atSite )
c              write (*,'( 10f10.4)') mag, rupLength, fltLength(iWidth), xl_F, P_SR_atSite

c             Find the number of rupture locations that rupture past the site
c             This is used to set the probability of the location (uniform)
              dx = 1.0  
              nLoc = int( (fltLength(iWidth) - rupLength)/ dx) 
              if ( nLoc .le. 1 ) then
                nLoc = 1
                rupLength = fltLength(iWidth)
              else
                dx = (fltLength(iWidth) - rupLength) / (nLoc-1)
              endif
                
              nLoc2 = 0
              nLoc2a = 0
              nloc2b = 0
              do iLoc=1,nLoc
c              set locations of the two ends of the rupture
               x1 = (iLoc-1)*dx
               x2 = x1 + rupLength
               if ( x1 .le. x_site0 .and. x2 .ge. x_site0) nLoc2 = nLoc2 + 1 
              enddo
c              write (80,'( 2i5)') nLoc2a, nLoc2b
              P_loc = 1./float(nLoc2)

c             Loop over rupture location (aleatory)
              do iLoc=1,nLoc

c              set locations of the two ends of the rupture
               x1 = (iLoc-1)*dx
               x2 = x1 + rupLength

c              skip if rupture is not past site            
                 if ( x1 .gt. x_site0 .or. x2 .lt. x_site0) goto 200

c              Set X/L for this rupture
               if ( xl_F .le. 0.5) then
                 X_over_L = (x_site0 - x1) / RupLength
               else
                 X_over_L = (x2 - x_site0) / RupLength
               endif

c              Fold the X_over_L to be in [0,0.5] range
               XL1 = X_over_L
               if ( XL1 .gt. 0.5 ) XL1 = 1.-XL1

c              Set XL bin
               iXLBin = int( (XL1-XL0)/XLBinStep ) + 1
       
c              Calculate median of Displacement (Dagg or Dp)
               call Calc_FDM (  jcalc_FDM(iFDM), Mag, iSOF, XL1, 
     1             muD, sigmaD, P_Gap, P_DP_zero, b2, mu_Dagg, sigma_Dagg )

c              Calculate sigma of Displacement (Dagg or Dp)
c              check for user-specified sigma
               if ( jcalc_sig(iFDM) .eq. 0 ) then 
                sigmaD   = fix_sigma(iFDM)
               else
                call Calc_FDM (  jcalc_Sig(iFDM),  Mag, iSOF, XL1, 
     1             muDX, sigmaD, P_GapX, P_DP_zeroX, b2, mu_Dagg, sigma_Dagg )
               endif
                 
c              add constant for epistemic unc in the median    
               muD = muD + deltamu(IFDM) 
                              
c              add constant for epitemic uncertainty in sigma
               sigmaD = sigmaD + delta_sig(iFDM)

c              compute rate * wt for for this scenario
               rate0 = rate(iMag,iBR,iWidth) * P_surf_rup * Prob_RL * P_SR_atSite 
     1             * P_loc * endBR_wt(iBR_total)
c               write (*,'( 10e12.3)') rate0, rate(iMag,iBR,iWidth) , P_surf_rup , Prob_RL , P_SR_atSite 
c     1             , P_loc , endBR_wt(iBR_total)
c               pause
               rate1 = rate(iMag,iBR,iWidth) * P_surf_rup * Prob_RL * P_SR_atSite 
     1             * P_loc 
               rate2 = rate(iMag,iBR,iWidth) * P_surf_rup * Prob_RL 
     1             * P_loc * endBR_wt(iBR_total)
               rate3 = rate(iMag,iBR,iWidth) * Prob_RL * P_SR_atSite 
     1             * P_loc * endBR_wt(iBR_total)
               rate4 = rate(iMag,iBR,iWidth) * P_surf_rup * Prob_RL * P_SR_atSite 
     1             * P_loc * endBR_wt4(iBR_total)
     
c              These sums are used to compute the mean PSR and P_rupatsite for each magnitude
               sum_b2(imag) = sum_b2(imag) + rate0
               sum_b3(imag) = sum_b3(imag) + rate2
               sum_b4(imag) = sum_b4(imag) + rate3
                
c              loop over each test value for hazard
               do iz=1,nZ

c               Set sigmatruncation level
                sigTrunc = 6.
      
c               compute the probability of exceedance of aggregate or primary rupture  
                ti = z(iz)**(0.3)
                Pz = pxceed3 ( muD, ti, sigmaD, sigTrunc ) 
                epsilonD = (ti - muD) / sigmaD
                     
c               set the epsilon bin (with                
                iEpsBin = int((epsilonD-Eps0) /EpsBinStep ) + 1
                if (iEpsBin .lt. 1 ) iEpsBin=1
                if (iEpsBin .gt. nEpsBin ) iEpsBin=nEpsBin
                
c               compute marginal hazard 
c               Add marginal hazard to the meanEps here
                marginalHaz = rate0 * Pz * (1.-P_Gap) * (1.- P_DP_zero) 
                marginalHaz_a = rate1 * Pz * (1.-P_Gap) * (1.- P_DP_zero) 
                marginalHaz_c = rate0 * Pz  * (1.- P_DP_zero) 
                marginalHaz_d = rate4 * Pz * (1.-P_Gap) * (1.- P_DP_zero) 
                meanEps(iz) = meanEps(iz) +  marginalHaz * epsilonD 
             
c               Add marginal hazard to hazard for a single mag (Mhaz) and to MeanMag and MeanXL           
                meanMag(iz) = meanMag(iz) +  marginalHaz * mag
                meanXL(iz) = meanXL(iz) +  marginalHaz * XL1
                MHaz(iMag,iz) = MHaz(iMag,iz) + marginalHaz
                deagg1(iz,iMagBin,iEpsBin) = deagg1(iz,iMagBin,iEpsBin) +  marginalHaz
                deagg2(iz,iMagBin,iXLBin) = deagg2(iz,iMagBin,iXLBin) +  marginalHaz

c               Add marginal hazard to the total mean hazard, haz_by_mag, and haz_by_branch
                Haz(iz) = Haz(iz) + marginalHaz
                Haz_mag(iz,iMag) = Haz_mag(iz,iMag) + marginalHaz
                Haz1(iz,iBR_total) = Haz1(iz,iBR_total) + marginalHaz_a
                Haz_FDM(iz,iFDM) = Haz_FDM(iz,iFDM) + marginalHaz_d
                
                if ( Haz(iz) .ge. 0. .and. Haz(iz) .le. 1. ) goto 300
                 write (*,'( 10e12.3)') marginalHaz, rate0, Pz, (1.-P_Gap), (1.- P_DP_zero) 
                 write (*,'( 10e12.3)') muD, ti, sigmaD, sigTrunc
                 stop 99
 300            continue

c               Add marginal hazard for hazard w/o gap              
                haz_noGap_prob(iZ) = haz_noGap_prob(iZ) + marginalHaz_c

c               set the epsilon bin for marginal disagg                
                iEpsBin1 = int((epsilonD-mineps1) /dEps1 ) + 1
                if (iEpsBin1 .lt. 1 ) iEpsBin1=1
                if (iEpsBin1 .gt. 61 ) iEpsBin1=61
                marg_Disag_eps(iz,iEpsBin1) = marg_Disag_eps(iz,iEpsBin1) + marginalHaz

c               set the Mag bin for marginal disagg                
                iMagBin1 = iMag
                marg_Disag_mag(iz,iMagBin1) = marg_Disag_mag(iz,iMagBin1) + marginalHaz

c               set the X/L bin for marginal disagg                
                iXLBin1 = int((XL1-minXL1) /dXL1 ) + 1
                if (iXLBin1 .lt. 1 ) iXLBin1=1
                if (iXLBin1 .gt. 61 ) iXLBin1=nXL1
                marg_Disag_XL(iz,iXLBin1) = marg_Disag_XL(iz,iXLBin1) + marginalHaz
          
               enddo
c              end loop over iz

c ------------------------- Distributed rupture hazard -------------

c              if site is on main fault (Rx=0) then skip distributed rupture calc
               if (Rx .eq. 0. ) goto 120

c              Set probability of distributed ruptures occurring along strike
c              within the along-strike dimension of the site (L_site)             
               P_distributed = 0.063 * L_site / 100.

c              Set probability that a distributed rupture is within the perpendicular to strike  
c              site dimension defined by the site width (
               call calc_PD_site ( w_site, Rx, iSOF, PD_site1, PD_site2)   

c              Loop over the aleatory variability of the aggregate rupture
c              (-3 to 3 epsilons)
               nAgg = 31 
               dEps = 0.2
               eps0D = -3. + dEps/2.
               do iAgg=1,nAgg

c               set epsilon for the D_agg
                epsD = eps0D + (iAgg-1) * dEps
                eps1D = epsD - dEps/2.
                eps2D = epsD + dEps/2.

c               compute d_agg (in M^0.3 units) and Probability for this epsilon                
                D_agg1 = mu_Dagg + eps*sigma_Dagg 
                Prob_D_agg = pxceed3 ( 0., eps1D, 1., 3. )- pxceed3 ( 0., eps2D, 1., 3. )

c               compute the total distributed displacement (m)
                if ( D_agg1 + b2 .gt. 0. ) then
                  D_Distrib_total = D_agg1**(1./0.3) - (D_agg1 + b2)**(1./0.3)
                else
                  D_Distrib_total = D_agg1**(1./0.3)
                endif

c               Loop over the number of distributed ruptures (Rd)
                nRd_max = 10
                sum33 = 0.
                do nRd = 1,nRd_max
                 call calc_N_distributed ( nRd, P_Rd, iSOF)

c                Random samples of the amplitudes
                 if ( nRd .eq. 1 ) then
                  nRatio_distrib = 1
                 else
                  nRatio_distrib = 20
                 endif

c                loop over number of amplitudes of distributed rupture for this
c                number of distributed ruptures
                 do iRatioD = 1, nRatio_distrib  
     
c                 Compute the distributed displacement ratio and the probability     
                  call calc_Amp_distrib ( nRd, iRatioD, DAMP_ratio, P_ratio_Damp)
                  D2 = (D_Distrib_total * DAMP_ratio)
c                  phi = D2_med*0.23
c                  tau = 0.201
c                  sigma_DT = sqrt(tau**2 + phi**2) 
                                  
c                 Compute the marginal hazard from distributed ruptures
c                 these are step functions
                  do iz=1,nZD
                   if ( D2 .gt. zd(iz) ) then
                     marginalHaz = rate0 * P_distributed * P_ratio_Damp * P_Rd * nRd
                     if ( iz .eq. 1 ) marghaz1 = marginalHaz
                   else
                     marginalHaz  = 0.
                   endif
                   HazS(1,iz) = HazS(1,iz) + marginalHaz * PD_site1
                   HazS(2,iz) = HazS(2,iz) + marginalHaz * PD_site2

                  enddo
c                 end loop over z values for distributed rupture
                  
                  write (*,'(4e12.3,i5,3e12.3 )') rate0, P_distributed, 
     1               P_ratio_Damp, P_Rd, nRd, margHaz1, D2, HazS(1,1)
                  sum33 = sum33 + P_distributed * P_ratio_Damp * P_Rd * nRd
                 enddo
c                end loop on the amp distributed                 
              
                enddo
c               end loop over number of distribued ruptures
                write (*,'( ''sum33:'',4e12.3)') sum33,pd_site1, pd_site2, HazS(1,1)
                

               enddo
c              end loop over D_agg variability

 120           continue
c              jump for FP=1 (zero amp for distributed)

c              skip to here if rupture is not past the site
 200           continue          

              enddo
c             end loop over rupture location

             enddo
c            end loop over rupture length variability        
    
            enddo
c           end loop over magnitude 
     
           enddo
c          end loops over logic tree branches for rate           

          enddo
c         end loop over Seismogenic thickness
       
         enddo    
c        end loop over SOF

        enddo
c       end loop over rupture length model

       enddo
c      end loop over PSR model    

      enddo
c     end loop over IFDM
      nEndBr = iBR_total
       
c ------------------------------------------------------------
      
c     Compute mean values for disaggregation    
      do iz=1,nZ
        if (  Haz(iz) .gt. 0. ) then
          meanMag(iz) = meanMag(iz) / Haz(iz)
          meanEps(iz) = meanEps(iz) / Haz(iz)
          meanXL(iz) = meanXL(iz) / Haz(iz)
        else
          meanMag(iz) = -99.
          meanEps(iz) =  -99.
          meanXL(iz) =  -99.
        endif
      enddo

c     write out disaggregation for M and eps bins  
      do iMagBin=1,nMagBin
       mag = mag0 + (iMagBin-1)*magBinStep
       
       do iEpsBin=1,nEpsBin
        eps = eps0 + (iEpsBin-1)*EpsBinStep
        do iz=1,nZ
         if (Haz(iz) .gt. 0. ) then
           deagg1(iz,iMagBin,iEpsBin) = deagg1(iz,iMagBin,iEpsBin) / Haz(iz)
         else
           deagg1(iz,iMagBin,iEpsBin) = 0.
         endif
        enddo
        write (61,'(4f8.2,100f8.4)') mag,mag+magBinStep, eps, eps+EpsBinStep,
     1          (deagg1(iz,iMagBin,iEpsBin),iz=1,nZ)
       enddo
      enddo

c     write out disaggregation for M, X/L bins      
      do iMagBin=1,nMagBin
       mag = mag0 + (iMagBin-1)*magBinStep

       do iXLBin=1,nXLBin
        XL1 = XL0 + (iXLBin-1)*XLBinStep
        do iz=1,nZ
         if (Haz(iz) .gt. 0. ) then
          deagg2(iz,iMagBin,iXLBin) = deagg2(iz,iMagBin,iXLBin) / Haz(iz)
         else
          deagg2(iz,iMagBin,iXLBin)  = 0.
         endif
        enddo
         write (61,'(4f8.2,100f8.4)') mag,mag+magBinStep, XL1, XL1+XLBinStep,
     1          (deagg2(iz,iMagBin,iXLBin),iz=1,nZ)
       enddo
      enddo
      
c     Write the marginal disagg for mag
      write (61,'( 2x,''  '')')
      write (61,'( 2x,''-----------------------------'')')
      write (61,'( 2x,''Marginal disagg for Mag'')')
      write (61,'( 18x,''mag'', 100f12.4)') (z(iz),iz=1,nz)
      
      do i=1,nMag1
        write (61,'( f10.2,100e12.3)') mag1(1,i), 
     1     (marg_disag_Mag(iz,i)/haz(iz),iz=1,nz)       
      enddo

c     Write the marginal disagg for X/L
      write (61,'( 2x,''  '')')
      write (61,'( 2x,''-----------------------------'')')
      write (61,'( 2x,''Marginal disagg for X/L'')')
      write (61,'( 18x,''mag'', 100f12.4)') (z(iz),iz=1,nz)
      do i=1,nXL1
        write (61,'( f10.2,100f10.4)') dXL1*(i-1), (marg_disag_XL(iz,i)/haz(iz),iz=1,nz)       
      enddo

c     Write the marginal disagg for eps
      write (61,'( 2x,''  '')')
      write (61,'( 2x,''-----------------------------'')')
      write (61,'( 2x,''Marginal disagg for Epsilon'')')
      write (61,'( 18x,''mag'', 100f12.4)') (z(iz),iz=1,nz)
      do i=1,nEps1
        write (61,'( f10.2,100f10.4)') minEps1 + dEps1*(i-1), (marg_disag_eps(iz,i)/haz(iz),iz=1,nz)       
      enddo

c     write the hazard by magnitude
      nMag2 = 0
      do iWidth=1,nWidth
        nMag2 = max(nMag2,nMag(iWidth))
      enddo
      do iMag=1,nMag2
        mag = mag1(1,iMag)
        write (60,'(''Marginal_Haz   '',2x,f5.2, 2x,100e12.3)') mag, 
     1            (mHaz(iMag,iz),iz=1,nz) 
      enddo
      
      write (*,'( 26x, 100f12.4)') (z(iz),iz=1,nz)
      write (*,'(''Mean_Total_Haz'',8x,100e12.3)') (Haz(iz),iz=1,nz) 
      write (*,'(''Mean_Mag      '',8x,100f12.2)') (meanMag(iz),iz=1,nz)
      write (*,'(''Mean_epsilon  '',8x,100f12.2)') (meanEps(iz),iz=1,nz)
      write (*,'(''Mean_X/L      '',8x,100f12.2)') (meanXL(iz),iz=1,nz)

c     Write mean hazard results
      write (60,'(''Mean_Total_Haz'',''     x    '',100e12.3)') (Haz(iz),iz=1,nz) 
      write (60,'(''Mean_Mag      '',''     x    '',100f12.2)') (meanMag(iz),iz=1,nz)
      write (60,'(''Mean_epsilon  '',''     x    '',100f12.2)') (meanEps(iz),iz=1,nz)
      write (60,'(''Mean_X/L      '',''     x    '',100f12.2)') (meanXL(iz),iz=1,nz)
            
c     write out the hazard without gap 
      write (60,'( /,''    z(m)     meanHaz     no_GapProb '')') 
      do iZ=1,nZ
        write (60,'(2x,f8.4, 5e12.3)') z(iz), Haz(iz), haz_noGap_prob(iZ)
      enddo    

c     write probability of rupture goes past the site and PSR by mag
      write (60,'( /,''----- Mean Prob_rup_past_site and PSR ----'')')
      write (60,'( '' Mag    P_rup_at_site  Mean_PSR  PSR_by_Model'')')
      do iMag=1,nMag2
        mag = mag1(1,iMag)
        if ( sum_b3(iMag) .gt. 0. ) then
          t1 = sum_b2(iMag)/sum_b3(iMag)
        else
          t1 = 0.
        endif
        if ( sum_b4(iMag) .gt. 0. ) then
          t2 = sum_b2(iMag)/sum_b4(iMag)
        else
          t2 = 0.
        endif
        write (60,'(f5.2, f10.4,4x,100f10.4)') mag, t1, t2,
     1       (PSR_save(iPSR,iMag),iPSR=1,nPSR)
      enddo
      
c     Write hazard by FDM and median disp given surface rup
      write (60,'( /,''----- Med disp given surface rupture at the site ----'')')
      write (60,'( '' iFDM   Med disp (m)'')')
      do iFDM=1,nFDM
       rate1 = haz_FDM(1,iFDM) / 2.
       z_med = 0.
       do iz=1,nz-1
         if ( haz_FDM(iz,iFDM) .ge. rate1 .and. haz_FDM(iz+1,iFDM) .le. rate1 ) then
           z_med = alog(z(iz)) + alog(z(iz+1)/z(iz)) * alog(rate1/haz_FDM(iz,iFDM))
     1           / alog(haz_FDM(iz+1,iFDM)/haz_FDM(iz,iFDM))
         endif
       enddo
       if ( z_med .ne. 0. ) then
         write (60,'( i5,f10.3)') iFDM, exp(z_med)
       else
         write (60,'( i5,2x,''add a larger z value'')') iFDM, exp(z_med)
       endif
      enddo
        
c     compute fractiles
      call Calc_frac ( haz1, endBR_wt, nEndBr, nZ, haz_frac )  

      write (30,'( '' ---------------------- Fractiles --------------------'')')
      write (30,'( 2x,''frac '',100f10.3)')  (z(iz),iz=1,nZ)
      do i=1,19
        write (30,'( 2x,f5.2,100e10.3)')  i*0.05, (haz_frac(iz,i),iz=1,nZ)
      enddo
      write (30,'( 2x,''Mean '',100e10.3)')  (haz_frac(iz,20),iz=1,nZ)
      
c ------------ end output for principal rupture hazard ---------------            
      
c --------- start output for distributed rupture hazard
      if ( rx .gt. 0. ) then

c       Loop over geologic information on surface faulting in the site region
        do ii =1,2
         write (*,'( /,'' ---------------------------------------------'')')
         write (60,'( /,'' ---------------------------------------------'')')
         if ( ii .eq. 1 ) write (60,'( 2x,''Distributed rupture hazard:'',
     1         '' Detailed geologic study found no faults in site region'')')
         if ( ii .eq. 2 ) write (60,'( 2x,''Distributed rupture hazard:'',
     2         '' Detailed geologic study in site region was not conducted'')')

         write (60,'( 100f12.3)') (z(iz),iz=1,nz)
         write (60,'(100e12.3))') (HazS(ii,iz),iz=1,nz) 
            
        enddo
      endif

 1000 write(*,'( 2x,''Normal termination'')')

      stop
      
      end
           



      subroutine test_Model
      real Mag, X_over_L, 
     1       mu_agg, mu_p, Dagg_sigma, DP_sigma,
     1       mu_agg_prime, mu_p_prime, Dagg_sigma_prime, DP_sigma_prime,
     3       P_Gap, P_DP_zero
      integer iSOF
      do i=1,10
        write (*,'( 2x,''Enter Mag, SOF, X/L'')')
        read (*,*) mag, iSOF, X_over_L

        call LA23_fDM ( Mag, iSOF, X_over_L, 
     1       mu_agg, mu_p, Dagg_sigma, DP_sigma,
     1       mu_agg_prime, mu_p_prime, Dagg_sigma_prime, DP_sigma_prime,
     3       P_Gap, P_DP_zero  )

        write (*,'( f10.3, 5x,''mu_agg          '', f10.3)') mu_agg**(10./3.)
        write (*,'( f10.3, 5x,''mu_P            '', f10.3)') mu_p**(10./3.)
        write (*,'( f10.3, 5x,''sigma_Dagg      '', f10.3)') dagg_sigma
        write (*,'( f10.3, 5x,''sigma_DP        '', f10.3)') dp_sigma
        write (*,'( f10.3, 5x,''mu_agg_prime    '', f10.3)') mu_agg_prime**(10./3.)
        write (*,'( f10.3, 5x,''mu_P_prime      '', f10.3)') mu_p_prime**(10./3.)
        write (*,'( f10.3, 5x,''sigma_Dagg_prime'', f10.3)') dagg_sigma_prime
        write (*,'( f10.3, 5x,''sigma_Dp_prime  '', f10.3)') dp_sigma_prime
        write (*,'( f10.3, 5x,''P(gap)          '', f10.3)') P_gap
        write (*,'( f10.3, 5x,''P(Dp=0)'', f10.3)') P_DP_zero
      enddo
      return
      end

c     -----------------------------

      subroutine CheckDim ( n1, nDIM, title )
      integer n1, nDIM
      character*20 title

      if ( n1 .gt. nDiM) then
        write (*,'( 2x,''Dimension error: Increase'',a20, ''to '',i5)') n1
        stop 99
      endif

      return
      end
