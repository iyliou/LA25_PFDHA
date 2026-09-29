
      subroutine rd_input ( nFDM, jcalc_FDM, jcalc_Sig, deltaMU, delta_sig, FDM_wt,
     1       nPSR, jcalc_PSR, PSR_wt, nRLM, jcalc_RLM, RLM_wt,
     2       nWidth, width_wt, nMag, mag1, nBR, rate_wt, rate,
     3       nSOF, iSOF1, SOF_wt, fltLength, seismoThick, dip, flt_width_DD,
     4       nZ, z, xl_F, Rx, w_site, L_site, F_map, iFlag_DD_prob, nZD, zD, 
     5       magBinStep,nMagBin, mag0, epsBinStep, nEpsBin, eps0,
     6       XLBinStep, nXLBin, XL0, fixed_sig, c1, c2, maxmag1, minmag1)
      implicit none
      include "rupHaz.inc"
      
      character*80 file1, fileout, dummy 
      real rate(MAX_NMAG,MAX_BR,MAXWIDTH), rate_wt(MAX_BR,MAXWIDTH)
      real mag1(MAXWIDTH,MAX_NMAG), maxmag1, minMag1
      real fltLength(MAXWIDTH), dip(MAXWIDTH), seismoThick(MAXWIDTH), flt_width_DD(MAXWIDTH)
      real PSR_wt(MAX_PSR), RLM_wt(MAX_RLM), SOF_wt(3)
      real z(MAXZ)
      real XL_F
      real width_wt(MAXWIDTH)
      real MagbinStep, EpsBinStep, mag0, eps0, XLBinStep, XL0
      real delta_sig0(MAX_FDM), delta_sig(MAX_FDM), delta_sig_wt(MAX_FDM)
      real deltaMU1(MAX_FDM), deltaMU2(MAX_FDM),deltaMU3(MAX_FDM), deltaMU(MAX_FDM)
      real FDM_wt(MAX_FDM), wt_mu1(MAX_FDM), wt_mu3(MAX_FDM)
      real fixed_sig(MAX_FDM), fixed_sig0(MAX_FDM)     
      real c1(MAX_PSR), c2(MAX_PSR)
      real sumwt

      integer jcalc_FDM(MAX_FDM), jcalc_Sig(MAX_FDM)
      integer jcalc_FDM0(MAX_FDM), jcalc_Sig0(MAX_FDM)
      integer i, iZ, nZ, k, i1, i2
      integer iMag, nMag(100)
      integer iSOF1(3), iSOF
      integer nSOF, nBr(MAXWIDTH), iBr
      integer nWidth, iWidth, nWidth1
      integer nMagBin, nEpsBin, nXLBin
      integer iPSR, iRLM, nPSR, nRLM
      integer jcalc_RLM(MAX_RLM), jcalc_PSR(MAX_PSR)
      integer iSig, nSig
      integer iFDM1, iFDm2, nFDM1, nFDM2, iFDM
      integer nFlt
      
c     distributed rupture param
      real zd(maxZ), Rx, W_site, L_site
      integer F_map, iFlag_DD_Prob
      integer  nFDM, nZd
      
      maxmag1 = 0.
      minmag1 = 10.

c     Read input file
      write (*,'( 2x,''Enter run file'')')
      read (*,'( a80)',err=1098) file1
      write (*,'( a80)') file1
      open (10,file=file1,status='old',err=1099)

c     open output files  
      read (10,'( a80)') dummy          
      read (10,'( a80)',err=2000) fileout
      open (60,file=fileout,status='unknown')
      read (10,'( a80)',err=2000) fileout
      open (30,file=fileout,status='unknown')
      read (10,'( a80)',err=2000) fileout
      open (61,file=fileout,status='unknown')      
      
c     read logic tree node for the FDM and constant terms
      read (10,'( a1)') dummy
      read (10,*,err=2001) nFDM1
      do iFDM1 = 1,nFDM1
        read (10,*,err=2002) jcalc_FDM0(iFDM1), deltaMU1(iFDM1), deltaMu2(iFDM1),
     1              jcalc_sig0(iFDM1), wt_mu1(iFDM1)
        fixed_sig0(iFDM1) = -99.
        write (*,'( 2x,''check jclac_sig'',i5)') jcalc_sig0(iFDM1)
        if ( jcalc_sig0(iFDM1) .le. 0 ) then
          write (*,'( 2x,''reading fixed sigma'', i5)') iFDM1
          backspace(10)
          read (10,*,err=2003) jcalc_FDM0(iFDM1), deltaMU1(iFDM1), deltaMu2(iFDM1),
     1              jcalc_sig0(iFDM1), wt_mu1(iFDM1), fixed_sig0(iFDM1)
          if ( fixed_sig0(iFDM1) .lt. 0. ) then
            write (*,'( 2x,''invalid fixed sigma'')')
            write (*,'( ''iFDM = '',i5,2x,''sigma='',f10.3)') iFDM1, fixed_sig0(iFDM1) 
          endif
          if ( fixed_sig0(iFDM1) .eq. 0. ) fixed_sig0(iFDM1) = 0.0001
        endif
        write (*,'( 2i5,f10.3)') iFDM1, jcalc_sig0(iFDM1), fixed_sig0(iFDM1)
      enddo
 
c     read logic tree node for epistemic uncertainty in non-ergodic terms
      read (10,'( a1)') dummy
      read (10,*,err=2004) nFDM2
      do iFDM2 = 1,nFDM2
        read (10,*,err=2005) deltaMU3(iFDM2), wt_mu3(iFDM2)
      enddo

c     Read logic tree node for adjustments to the sigma for the FDM
      read (10,'( a1)') dummy
      read (10,*,err=2006) nSig
      do iSig = 1,nSig
        read (10,*,err=2007) delta_Sig0(iSig), delta_sig_wt(iSig)
      enddo
      
c     combine logic trees for the FDM into one array
      i = 1
      sumwt = 0.
      do iFDM1 = 1,nFDM1
       do iFDM2 = 1,nFDM2
        do iSig=1,nSig
         deltaMu(i) = deltaMu1(iFDM1) + deltaMu2(iFDM1) + deltaMu3(iFDM2)
         FDM_wt(i) = wt_mu1(iFDM1) * wt_mu3(iFDM2) * delta_sig_wt(iSig)
         jcalc_FDM(i) = jcalc_FDM0(iFDM1)
         jcalc_sig(i) = jcalc_sig0(iFDM1)
         fixed_sig(i) = fixed_sig0(iFDM1)
         delta_sig(i) = delta_Sig0(iSig)
         sumwt = sumwt + FDM_wt(i)
         i = i + 1
        enddo
       enddo
      enddo
      nFDM = i-1
      write (*,'( 2x,''check sum FDM_wt'',f10.4)') sumwt      

c     Set the probability for surface rupture model
      read (10,'( a1)') dummy
      sumwt = 0.
      read (10,*,err=2008) nPSR
      do iPSR=1,nPSR
        read (10,*,err=2009) jcalc_PSR(iPSR), PSR_wt(iPSR)
        if ( jcalc_PSR(iPSR) .eq. 1 ) then
          backspace (10)
          read (10,*,err=2009) jcalc_PSR(iPSR), PSR_wt(iPSR), c1(iPSR), c2(iPSR)
        endif
        sumwt = sumwt + PSR_wt(iPSR)
      enddo
      write (*,'( 2x,''check sum PSR_wt'',f10.4)') sumwt

c     Set the surface rupture length model
      read (10,'( a1)') dummy
      sumwt = 0.
      read (10,*,err=2010) nRLM
      do iRLM=1,nRLM
        read (10,*,err=2011) jcalc_RLM(iRLM), RLM_wt(iRLM)
        sumwt = sumwt + RLM_wt(iRLM)
      enddo
      write (*,'( 2x,''check sum RLM_wt'',f10.4)') sumwt
           
c     open out2 file from Haz45
      read (10,'( a1)') dummy
      read (10,'( a80)',err=2012) fileout
      write (30,'( a80)') fileout
      write (30,'( 2x,''out2 file name'')')
      open (20,file=fileout,status='old',err=2013)
      
c     Read out2 file      
      do i=1,2
        read (20,'( a80)') dummy
      enddo
      read (20,*) nFlt
      if ( nFlt .ne.1 ) then
        write (*,'( 2x,''Current version only setup for 1 flt'')')
        stop 99
      endif
      read (20,*,err=2014) nWidth
      read (20,*,err=2015) (width_wt(k),k=1,nWidth) 
c      write (*,'( i5,e12.3)')     (k, width_wt(k),k=1,nWidth) 
      
      do iWidth=1,nWidth
        read (20,'( a80)') dummy
        read (20,*,err=2016) nMag(iWidth), nBR(iWidth)
        read (20,*,err=2017) (rate_wt(iBR,iWidth),iBR=1,nBR(iWidth)) 
c        write (*,'(10e12.3)') (rate_wt(iBR,iWidth),iBR=1,nBR(iWidth)) 
        
        sumwt = 0.
        do iBR=1,nBR(iWidth)
          sumwt = sumwt + rate_wt(iBR,iWidth)
        enddo
        write (*,'( 2x,''check sum rate_wt for iWidth'',i5,f10.4)') iWidth, sumwt
c        write (*,'( i5)') nBR(iWidth)

        do iMag=1,nMag(iWidth)
          read (20,*,err=2018) i1, i2, mag1(iWidth,iMag)
     1       ,(rate(iMag,iBR,iWidth),iBR=1,nBR(iWidth))
          maxmag1 = max( maxmag1,mag1(iWidth,iMag))
          minmag1 = min( minmag1,mag1(iWidth,iMag))
c          write (*,'( 2i5,2f10.3,e12.3)') iMag,nMag(iWidth), maxmag1, minmag1, rate(iMag,1,iWidth)
        enddo
      enddo
            
c     read SOF logic tree  (SOF def:  -1 NML, 0 SS, 1 RV)
      read (10,'( a80)') dummy
      sumwt = 0.
      read (10,*,err=2019) nSOF
      do iSOF=1,nSOF
        read (10,*,err=2020) iSOF1(iSOF),SOF_wt(iSOF)
        sumwt = sumwt + SOF_wt(iSOF)
      enddo
      write (*,'( 2x,''check sum SOF_wt'',f10.4)') sumwt
      
c     Read flt geometry
      read (10,'( a80)') dummy
      read (10,*,err=2021) nWidth1
      if ( nWidth1 .ne. nWidth) then
        write (*,'( 2x,''inp0ut error, Check nWight from out2 file'')')
        stop 99
      endif
      do iWidth=1,nWidth
        read (10,*,err=2022) fltLength(iWidth), seismoThick(iWidth), dip(iWidth), width_wt(iWidth)
        flt_Width_DD(iWidth) = seismoThick(iWidth) / sin(dip(iWidth)*3.1415926/180.) 
      enddo
      write (*,'( 2x,''end of read flt geometry'')') 
      

c     Read location along rupture and dist perpendicular to strike   
      read (10,'( a80)') dummy                  
      read (10,*,err=2025) xl_F, Rx

c     If Rx ne 0, then read site dimension and z values for distributed hazard
      if ( Rx .ne. 0. ) then
        read (10,*) w_site, L_site
c        read (10,*) iFlag_DD_prob
        read (10,'( a80)') dummy      
        read (10,*) nZD 
        read (10,*) (zD(iz),iz=1,nzD)
      endif
      write (*,'( 2x,''end of read site'')') 

c     Read z values for principal disp hazard
      read (10,'( a80)') dummy          
      read (10,*,err=2023) nZ
      read (10,*,err=2024) (z(iz),iz=1,nz)

c     read Disaggregation inputs
      read (10,'( a80)') dummy      
      read (10,*,ERR=2026) magBinStep,nMagBin, mag0
      read (10,*,ERR=2027) epsBinStep, nEpsBin, eps0
      read (10,*,ERR=2028) XLBinStep, nXLBin, XL0
      if ( nMagBin .gt. MAXMagBin ) then
       write (*,'( 2x,''increase dimension for MaxMagBin'')')
       stop 99
      elseif ( nEpsBin .gt. MAXEpsBin ) then
       write (*,'( 2x,''increase dimension for MaxEpsBin'')')
       stop 99
      elseif ( nXLBin .gt. MaxXLBin ) then
       write (*,'( 2x,''increase dimension for MaxXLBin'')')
       stop 99
      endif
      
c     Echo inputs to output file 
    
      write (60,'( 2x,''----FDM----'')')
      write (60,'( ''    Med_index Sig_index deltaMu1  deltaMu2  FDM_wt   Fix_sigma'')')
      do iFDM1=1,nFDM1
        write (60,'(2i10,10f10.4)') jcalc_FDM0(iFDM1),  jcalc_Sig0(iFDM1), 
     1      deltaMU1(iFDM1), deltaMU2(iFDM1), 
     1       wt_mu1(iFDM1) , fixed_sig0(iFDM1)
      enddo

      write (60,'( /,2x,''----delta_NE_Epistemic----'')')
      write (60,'( ''    deltaMu3     wt   '')')
      do iFDM2=1,nFDM2
        write (60,'(210f10.4)') deltaMU3(iFDM2), wt_mu3(iFDM2)
      enddo

      write (60,'( /,2x,''----sigma_Epistemic----'')')
      write (60,'( ''    delta_sig     wt   '')')
      do iSig=1,nSig
        write (60,'(2f10.4)') delta_Sig0(iSig), delta_sig_wt(iSig)
      enddo

      write (60,'( /,2x,''----Combined FDM----'')')
      write (60,'( ''    Med_index Sig_index Net_deltaMu  delta_sig  net_wt   Fix_sigma'')')
      do iFDM=1,nFDM
        write (60,'(2i10,3x,10f10.4)') jcalc_FDM(iFDM),  jcalc_Sig(iFDM), 
     1      deltaMU(iFDM), delta_sig(iFDM), 
     1      FDM_wt(iFDM), fixed_sig(iFDM)
      enddo


      write (60,'( /,2x,''---- PSR Model----'')')
      write (60,'( ''     PSR_index   wt'')')
      do iPSR=1,nPSR
        write (60,'(i10,f10.3)') jcalc_PSR(iPSR), PSR_wt(iPSR)
      enddo
      
      write (60,'( /,2x,''---- RL Model ---'')')
      write (60,'( ''     RLM_index   wt'')')
      do iRLM=1,nRLM
        write (60,'(i10,f10.3)') jcalc_RLM(iRLM), RLM_wt(iRLM)        
      enddo

      write (60,'( /,2x,''---- width wts ----'')')
      write (60,'(''      iWidth   wt '' )') 
      do iWidth=1,nWidth
        write (60,'( i10,f10.4)') iWidth, width_wt(iWidth)
      enddo
      
      write (60,'( /,2x,''---- rate by width ----'')')
      write (60,'(''  iWidth   IMAG     Mag    rates_by_branch '' )') 
      do iWidth=1,nWidth
        do iMag=1,nMag(iWidth)
          write (60,'( 2i7, f10.2,100e12.3)') iWidth, IMAG, mag1(iWidth,imag),
     1          (rate(iMag,iBR,iWidth),IBR=1,nBR(iWidth))
        enddo
      enddo

      write (60,'( /,2x,''---- SOF ----'')')
      write (60,'( ''        SOF      wt'')')
      do iSOF=1,nSOF
        write (60,'( i10,f10.3)') iSOF1(iSOF), SOF_wt(ISOF)
      enddo
      
      write (60,'( /,2x,''---- flt geometry by width ----'')')
      write (60,'( ''      iWidth  Flt_Len  SeismoThick   Dip      wt'')')
      do iWidth=1,nWidth
        write (60,'( i10,4f10.2)') iWidth, fltLength(iWidth), 
     1       seismoThick(iWidth), dip(iWidth), width_wt(iWidth)
      enddo
   
      write (60,'( /,''----  site location ---- '')')
      write (60,'( ''     XF/LF      Rx'')')
      write (60,'( 2f10.3)')  xl_F, Rx
                             
      
      return
      
 1098 write (*,'( 2x,''error in run file name'')') 
      goto 2500
 1099  write (*,'( 2x,''run file not found'')') 
      goto 2500
     
 2000 write (*,'( 2x,''error in input: output file name'')') 
      goto 2500
 2001 write (*,'( 2x,''error in input: nFDM'')') 
      goto 2500
 2002 write (*,'( 2x,''error in input: FDM model'')') 
      goto 2500
 2003 write (*,'( 2x,''error in input: fixed sigma'')') 
      goto 2500
 2004 write (*,'( 2x,''error in input: nDelta_NE'')') 
      goto 2500
 2005 write (*,'( 2x,''error in input: delta_NE epistemic'')') 
      goto 2500
 2006 write (*,'( 2x,''error in input: nSigma epi'')') 
      goto 2500
 2007 write (*,'( 2x,''error in input: sigma epi'')') 
      goto 2500
 2008 write (*,'( 2x,''error in input: nPSR models'')') 
      goto 2500
 2009 write (*,'( 2x,''error in input: PSR models, wts'')') 
      goto 2500
 2010 write (*,'( 2x,''error in input: nRL models'')') 
      goto 2500
 2011 write (*,'( 2x,''error in input: RL models, wts'')') 
      goto 2500
 2012 write (*,'( 2x,''error in out2 file name'')') 
      goto 2500
 2013 write (*,'( 2x,''error opening out2 file'')') 
      goto 2500
 2014 write (*,'( 2x,''error in out2 file: nWidth'')') 
      goto 2501
 2015 write (*,'( 2x,''error in out2 file: wts for Width'')') 
      goto 2501
 2016 write (*,'( 2x,''error in out2 file: nMag'')') 
      goto 2501
 2017 write (*,'( 2x,''error in out2 file: wts for branches'')') 
      goto 2501
 2018 write (*,'( 2x,''error in out2 file: mag, rates'')') 
      goto 2501
 2019 write (*,'( 2x,''error in input: nSOF'')') 
      goto 2500
 2020 write (*,'( 2x,''error in input: SOF and wts'')') 
      goto 2500
 2021 write (*,'( 2x,''error in input: nWidth '')') 
      goto 2500
 2022 write (*,'( 2x,''error in input: fault Len, seismoThick, dip '')') 
      goto 2500
 2023 write (*,'( 2x,''error in input: nDisp for haz'')') 
      goto 2500
 2024 write (*,'( 2x,''error in input: disp for haz'')') 
      goto 2500
 2025 write (*,'( 2x,''error in input: site XF/LF, Rx '')') 
      goto 2500
 2026 write (*,'( 2x,''error in input: Deagg bins for Mag'')') 
      goto 2500
 2027 write (*,'( 2x,''error in input: Deagg bins for epsilon'')') 
      goto 2500
 2028 write (*,'( 2x,''error in input: Deagg bins for X/L'')') 
      goto 2500
      
 2500 backspace (10)
      read (10,'( a80)') dummy
      write (*,'( 2x,1x,a80)') dummy
      stop

 2501 backspace (20)
      read (20,'( a80)') dummy
      write (*,'( 2x,1x,a80)') dummy
      stop
      
      end

 
