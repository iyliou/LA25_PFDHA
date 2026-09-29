      subroutine nonerg_sigma ( iSigModel, mag, XL, disp, iSOF, sigma, iSig_form )
      
      integer iSigModel, iSig_form, iSOF
      real mag, XL, sigma, disp
      
c     Liou and ABrahamson (2025) non-ergodic model 1
c     for magpdf sigma = 0.10
      if ( iSigModel .eq. 101 ) then 
        call LA25_NE_sigma1( disp, sigma  )
        return
      endif
      
c     Liou and ABrahamson (2025) non-ergodic model 2
c     for magpdf sigma = 0.15
      if ( iSigModel .eq. 102 ) then 
        call LA25_NE_sigma2( disp, sigma  )
        return
      endif

      write (*,'( 2x,''bad Nonergodic sigma model'',i5)') iSigModel
      stop 99
      end
      
c    -----------------------------

      subroutine LA25_NE_sigma1( DS_bar, sigma  )
      
      real DS_bar, sigma_MX
      real sigma_MX_NE, sigma_NE
      real d1, a0 a1

c     DS_bar is the arithmetic mean of the displacements at the site

c     set coeff for sigma_MX for sigmaM=0.1 (Table 5)
      d1 = 0.38
      a0 = 0.087
      a1 = -0.0178
      
c     Compute sigma_MX_NE (eq 22)
      if ( DS_bar .lt. 1.44 ) then
        sigma_MX_NE = 0.17
      else
        sigma_MX_NE = 0.17 + 0.0325 * ( alog(DS_bar/1.44))**2
      endif

c     compute sigma_MX (eq 24)      
      if ( DS_bar .lt. d1 ) then
        sigma_MX = a0
      else
        sigma_MX = a0 + a1 * (alog(DS_bar/d1))**2
      endif

c     Compute nonergodic sigms (eq 20)
      sigma_NE = sqrt( sigma_MX_NE**2 - sigma_MX**2)
      return
      end

c    -----------------------------

      subroutine LA25_NE_sigma2( DS_bar, sigma  )
      
      real DS_bar, sigma_MX
      real sigma_MX_NE, sigma_NE
      real d1, a0 a1

c     DS_bar is the arithmetic mean of the displacements at the site

c     set coeff for sigma_MX for sigmaM=0.15 (Table 5)
      d1 = 0.383
      a0 = 0.0994
      a1 = -0.0206
      
c     Compute sigma_MX_NE (eq 22)
      if ( DS_bar .lt. 1.44 ) then
        sigma_MX_NE = 0.17
      else
        sigma_MX_NE = 0.17 + 0.0325 * ( alog(DS_bar/1.44))**2
      endif

c     compute sigma_MX (eq 24)      
      if ( DS_bar .lt. d1 ) then
        sigma_MX = a0
      else
        sigma_MX = a0 + a1 * (alog(DS_bar/d1))**2
      endif

c     Compute nonergodic sigms (eq 20)
      sigma_NE = sqrt( sigma_MX_NE**2 - sigma_MX**2)
      return
      end

