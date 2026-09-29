      subroutine Set_seg_location ( nSeg, Li1, LT_flt, LT_Ratio, xSeg, nCase_SegLoc)
      implicit none
      integer nSeg, nCase_SegLoc
      integer nx, nx1, nx2, ix, ix1, ix2, jx, n, i
      real Li1(10), LT_flt, LT_Ratio, xSeg(4,10000)
      real x1_seg1, x2_seg1, x1_seg2, x2_seg2, x1_seg3, x2_seg3, x1_seg4, x2_seg4
      real dx, xx

      dx = LT_flt * 0.02
c      if ( dx .gt. 1. ) dx = 1.

c     Simple cases for nSeg=1 or nSeg=2      
      if ( nSeg .eq. 1 ) then
        xSeg(1,1) = 0.
        nCase_segLoc = 1
      elseif ( nSeg .eq. 2 ) then
        xSeg(1,1) = 0.0
        xSeg(2,1) = LT_Flt - Li1(2)
        nCase_segLoc = 1
      endif
      
c     Case for nSeg=3 and LT_ratio < 1  (gaps)
      if ( nSeg .eq. 3 .and. LT_ratio .le. 1. ) then
        x2_seg1 = Li1(1) 
        x1_seg3 = LT_flt - Li1(3)

        nx = int((x1_seg3 - x2_seg1 - Li1(2))/dx)
c        write (*,'( 5f10.3,i5)') x1_seg3, x2_seg1 , Li1(2), dx, ((x1_seg3 - x2_seg1 - Li1(2))/dx), nx
        if ( nx .lt. 1 ) nx = 1
        
        do ix=1,nx
          xSeg(1,ix) = 0.
          xSeg(3,ix) = LT_flt - Li1(3)
          xSeg(2,ix) = x2_seg1 + dx*(ix-1)
        enddo
        if ( xSeg(2,nx)+Li1(2) .gt. x1_seg3 ) xSeg(2,nx) = x1_seg3-Li1(2)
        nCase_segLoc = nx
c        write (*,'( 3f10.3,i5)') Li1(1), Li1(2), li1(3), nx
c        pause
      endif

c     Case for nSeg=3 and LT_ratio > 1  (overlap)
      if ( nSeg .eq. 3 .and. LT_ratio .gt. 1. ) then
        x1_seg1 = 0.
        x2_seg1 = Li1(1)
        x1_seg3 = LT_flt - Li1(3)
        x2_seg3 = LT_flt

        nx1 = int((LT_flt - Li1(2))/dx)
        if ( nx1 .lt. 1 ) nx1 = 1
        nx2 = int((LT_flt - Li1(3))/dx)
        if ( nx2 .lt. 1 ) nx1 = 2
        
        jx = 0
        do ix1=1,nx1
         x1_seg2 = 0. + dx*(ix1-1)
         x2_seg2 = x1_seg2 + Li1(2)

c        check that there are no gap, if so, keep this case
         n = int(LT_flt/dx)
         do i=1,n
           xx = (i-1)*dx
c           write (*,'( 7f10.2)') xx, x1_seg1, x2_seg1, x1_seg2, x2_seg2, x1_seg3, x2_seg3
           if ( x1_seg1 .le. xx .and. x2_seg1 .ge. xx ) goto 200
           if ( x1_seg2 .le. xx .and. x2_seg2 .ge. xx ) goto 200
           if ( x1_seg3 .le. xx .and. x2_seg3 .ge. xx ) goto 200
c           write (*,'( ''gap'',7f10.2)') xx, x1_seg1, x2_seg1, x1_seg2, x2_seg2, x1_seg3, x2_seg3
c            there is a gap, so skip this case
             goto 210
 200       continue
         enddo

c        keep this case
         jx = jx + 1
         if ( jx. gt. 10000) then
           write (*,'( 2x,''maxcase exceeded'',i6)') jx
           stop 99
         endif
         xSeg(1,jx) = 0.
         xSeg(2,jx) = x1_seg2
         xSeg(3,jx) = x1_seg3
 210     continue
        enddo
        nCase_segLoc = jx
        if (nCase_segLoc .eq. 0 ) then
          jx = 1 
          xSeg(1,jx) = 0.
          xSeg(2,jx) = Li1(1)
          xSeg(3,jx) = LT_flt - Li1(3)
          nCase_segLoc = 1
        endif
      endif

c     Case for nSeg=4 and LT_ratio < 1  (gaps)
      if ( nSeg .eq. 4 .and. LT_ratio .le. 1. ) then
        x1_seg1 = 0.
        x2_seg1 = Li1(1)
        x1_seg4 = LT_flt - Li1(4)
        x2_seg4 = LT_flt

        nx1 = int((x1_seg4-x2_seg1-Li1(3)-Li1(2))/dx)
        if ( nx1 .lt. 1 ) nx1 = 1
        
        jx = 0
        do ix1=1,nx1
         x1_seg2 = x2_seg1 + dx*(ix1-1)
         x2_seg2 = x1_seg2 + Li1(2)
         nx2 = int((x1_seg4-x2_seg2-Li1(3))/dx)
        if ( nx2 .lt. 1 ) nx2 = 1
         do ix2=1,nx2
           x1_seg3 = x2_seg2 + dx*(ix2-1)
           jx = jx + 1
           xSeg(1,jx) = x1_seg1 
           xSeg(4,jx) = x1_seg4
           xSeg(2,jx) = x1_seg2
           xSeg(3,jx) = x1_seg3 
         enddo
        enddo
c        if ( xSeg(3,jx)+Li1(3) .gt. x1_seg4 ) xSeg(3,jx) = x1_seg4-Li1(3)
        nCase_segLoc = jx
      endif

c     Case for nSeg=4 and LT_ratio > 1  (overlap)
      if ( nSeg .eq. 4 .and. LT_ratio .gt. 1. ) then
        x1_seg1 = 0.
        x2_seg1 = Li1(1)
        x1_seg4 = LT_flt - Li1(4)
        x2_seg4 = LT_flt
        
        nx1 = int((LT_flt - Li1(2))/dx)
        if ( nx1 .lt. 1 ) nx1 = 1
        nx2 = int((LT_flt - Li1(3))/dx)
        if ( nx2 .lt. 1 ) nx1 = 2
        
        jx = 0
        do ix1=1,nx1
         x1_seg2 = 0. + dx*(ix1-1)
         x2_seg2 = x1_seg2 + Li1(2)
         do ix2=1,nx2
           x1_seg3 = 0. + dx*(ix1-1)
           x2_seg3 = x1_seg3 + Li1(3)

c          check that there are no gap, if so, keep this case
           n = int(LT_flt/dx)
           do i=1,n
             xx = (i-1)*dx
             if ( x1_seg1 .le. xx .and. x2_seg1 .ge. xx ) goto 100
             if ( x1_seg2 .le. xx .and. x2_seg2 .ge. xx ) goto 100
             if ( x1_seg3 .le. xx .and. x2_seg3 .ge. xx ) goto 100
             if ( x1_seg4 .le. xx .and. x2_seg4 .ge. xx ) goto 100

c            there is a gap, so skip this case
c             write (*,'( ''gap'',9f10.2)') xx, x1_seg1, x2_seg1, x1_seg2, x2_seg2,
c     1        x1_seg3, x2_seg3,  x1_seg4, x2_seg4
             goto 110
 100         continue
           enddo
           
           jx = jx + 1
           if ( jx. gt. 10000) then
             write (*,'( 2x,''maxcase exceeded'',i6)') jx
             stop 99
           endif
           xSeg(1,jx) = 0.
           xSeg(2,jx) = x1_seg2
           xSeg(3,jx) = x1_seg3
           xSeg(4,jx) = x1_seg4
 110       continue
         enddo
        enddo
        nCase_segLoc = jx
        
        if (nCase_segLoc .eq. 0 ) then
          jx = 1 
          xSeg(1,jx) = 0.
          xSeg(2,jx) = Li1(1)
          xSeg(3,jx) = Li1(1) + Li1(2)
          xSeg(4,jx) = LT_flt-Li1(4)
          nCase_segLoc = 1
        endif
      endif
      
      return
      end

c --------------------

      subroutine calc_N_distributed ( iRd, P_Rd, SOF)
      implicit none
      integer iRd, SOF
      real P_Rd
      real SS(15), RV(15), NML(15)

      data NML / 0.547, 0.2478, 0.1122, 0.0508, 0.023, 0.0104, 0.0047,
     1          0.0021, 0.001, 0.0004, 0.0002, 0.0001, 0, 0, 0 /
      data SS / 0.304, 0.2116, 0.1473, 0.1025, 0.0713, 0.0497, 0.0346, 
     1 0.0241, 0.0167, 0.0117, 0.0081, 0.0056, 0.0039, 0.0027, 0.0019 /
      data RV / 0.809, 0.1545, 0.0295, 0.0056, 0.0011, 0.0002, 0, 0, 0,
     1           0, 0, 0, 0, 0, 0 /

      if (SOF .eq. 0 ) P_Rd = SS(iRd)
      if (SOF .eq. 1 ) P_Rd = RV(iRd)
      if (SOF .eq. -1 ) P_Rd = NML(iRd)
      
      return
      end    

c -------------------------

      subroutine calc_Amp_distrib ( nRd, iRatio, ratio, P_ratio)
      implicit none
      integer iRatio, nRd
      real ratio, P_ratio
      real ratio1(20), P_NRD2(20), P_NRD3(20), P_NRD4(20), P_NRD5(20)

      data ratio1 / 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325,
     1              0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 
     2              0.725, 0.775, 0.825, 0.875, 0.925, 0.975 /
      data P_NRD2 / 0.0058, 0.0135, 0.0198, 0.0253, 0.0304, 0.0351, 0.0396,
     1              0.0439, 0.048, 0.0519, 0.0556, 0.0592, 0.0626, 0.0658, 
     2              0.0688, 0.0716, 0.0741, 0.0762, 0.0774, 0.0755 /
      data P_NRD3 / 0.0182, 0.0348, 0.0449, 0.0523, 0.0578, 0.0618, 0.0647, 
     1              0.0665, 0.0673, 0.0672, 0.0662, 0.0645, 0.0618, 0.0583, 
     2              0.0539, 0.0485, 0.042, 0.0341, 0.0243, 0.011 /
      data P_NRD4 / 0.0488, 0.0746, 0.0842, 0.0878, 0.0878, 0.0853, 0.0812,
     1              0.0757, 0.0694, 0.0624, 0.0551, 0.0475, 0.0399, 0.0325, 
     2              0.0253, 0.0186, 0.0126, 0.0074, 0.0033, 0.0006 /
      data P_NRD5 / 0.0903, 0.1252, 0.1283, 0.1211, 0.1092, 0.0951, 0.0805,
     1              0.0663, 0.0531, 0.0412, 0.0309, 0.0223, 0.0154, 0.0099, 
     2              0.0059, 0.0032, 0.0014, 0.0005, 0.0001, 0 /     

      if (nRd .eq. 1 ) then
        ratio = 1.0
        P_ratio = 1.0
      else 
        ratio = ratio1(iRatio)
        if ( nRd .eq. 2 ) P_ratio = P_NRD2(iRatio)
        if ( nRd .eq. 3 ) P_ratio = P_NRD3(iRatio)
        if ( nRd .eq. 4 ) P_ratio = P_NRD4(iRatio)
        if ( nRd .ge. 5 ) P_ratio = P_NRD5(iRatio)
      endif

c     Scale for nRd > 5 (to get the mean of the sum to be unity)
      if ( nRd .gt. 5 ) then
        ratio = ratio * nRd / 5.
      endif
      
c      write (*,'( ''p_ratio'',e12.3)') P_ratio
      
      return
      end    

c --------------------

      subroutine calc_PD_site ( site_width, Rx, SOF, PD_site1, PD_site2)   
      implicit none
      real site_width, Rx, PD_site1, PD_site2 
      integer SOF
      real pi, pdf
      real sigma_SS, sigma_NML, sigma_RV, phi_SS, phi_NML, phi_RV
      real sigma, phi
      real Rmax, pdf1, pdf2, part1, part2
            
c     set sigma (m) for the normal distribution part      
      sigma_SS = 1515.
      sigma_NML = 761.
      sigma_RV = 647.

c     set parameters for the uniform part
      phi_NML = 0.024
      phi_SS = 0.258
      phi_RV = 0.058
      Rmax = 20000.

      if ( SOF .eq. -1 ) then
        sigma = sigma_NML
        phi = phi_NML
      elseif ( SOF .eq. 0 ) then
        sigma = sigma_SS
        phi = phi_SS
      elseif ( SOF .eq. 1 ) then
        sigma = sigma_RV
        phi = phi_RV
      else
        write (*,'( 2x,''bad SOF in calc_PD_site'', i5)') SOF
        stop 99
      endif
      
c     first term (half normal)
      pi = 3.1415926
      pdf1 = exp(-(Rx**2)/(2.*sigma**2)) / (sqrt(2.*pi)* sigma)
      
c     set second term (uniform)
      if ( Rx .lt. Rmax ) then
        pdf2 = 1./ Rmax
      else
        pdf2 = 0.
      endif

c     Set PD_site if geologic studies show no surface faulting in the site region
      PD_site1 = pdf1 * (1-phi) * site_width

c     Set PD_site if there are no geologic studies of the surface faulting in the site region
c     set amplitudes of the uniform and normal parts to have a total pdf that sums to unity
      part2 = pdf2 * phi
      part1 = pdf1 * (1-phi)
      PD_site2 = (part1 + part2) * site_width
      
c      write (*,'( 20e12.3)') sigma, phi, Rx, pdf1, pdf2, phi, part1, part2, PD_site1, PD_site2
c      pause

      return
      end 
c ------------------------------
      subroutine setOrder (nSeg, nOrder, order )
      integer order(100,5), nSeg, nORder

      if ( nSeg .eq. 1 ) then
          nOrder = 1
          order(1,1) = 1
      elseif ( nSeg .eq. 2 ) then
          nOrder = 2
          order(1,1) = 1
          order(2,1) = 2
          order(1,2) = 2
          order(2,2) = 1
      elseif ( nSeg .eq. 3 ) then
          nOrder = 6
          order(1,1) = 1
          order(1,2) = 2
          order(1,3) = 3
          order(2,1) = 1
          order(2,2) = 3
          order(2,3) = 2
          order(3,1) = 2
          order(3,2) = 1
          order(3,3) = 3
          order(4,1) = 2
          order(4,2) = 3
          order(4,3) = 1
          order(5,1) = 3
          order(5,2) = 1
          order(5,3) = 2
          order(6,1) = 3
          order(6,2) = 2
          order(6,3) = 1
      else
          nOrder = 24
          order(1,1) = 1
          order(1,2) = 2
          order(1,3) = 3
          order(1,4) = 4
          order(2,1) = 1
          order(2,2) = 2
          order(2,3) = 4
          order(2,4) = 3
          order(3,1) = 1
          order(3,2) = 3
          order(3,3) = 2
          order(3,4) = 4
          order(4,1) = 1
          order(4,2) = 3
          order(4,3) = 4
          order(4,4) = 2
          order(5,1) = 1
          order(5,2) = 4
          order(5,3) = 2
          order(5,4) = 3
          order(6,1) = 1
          order(6,2) = 4
          order(6,3) = 3
          order(6,4) = 2

          order(7,1) = 2
          order(7,2) = 1
          order(7,3) = 3
          order(7,4) = 4
          order(8,1) = 2
          order(8,2) = 1
          order(8,3) = 4
          order(8,4) = 3
          order(9,1) = 2
          order(9,2) = 3
          order(9,3) = 1
          order(9,4) = 4
          order(10,1) = 2
          order(10,2) = 3
          order(10,3) = 4
          order(10,4) = 1
          order(11,1) = 2
          order(11,2) = 4
          order(11,3) = 1
          order(11,4) = 3
          order(12,1) = 2
          order(12,2) = 4
          order(12,3) = 3
          order(12,4) = 1

          order(13,1) = 3
          order(13,2) = 1
          order(13,3) = 2
          order(13,4) = 4
          order(14,1) = 3
          order(14,2) = 1
          order(14,3) = 4
          order(14,4) = 2
          order(15,1) = 3
          order(15,2) = 2
          order(15,3) = 1
          order(15,4) = 4
          order(16,1) = 3
          order(16,2) = 2
          order(16,3) = 4
          order(16,4) = 1
          order(17,1) = 3
          order(17,2) = 4
          order(17,3) = 1
          order(17,4) = 2
          order(18,1) = 3
          order(18,2) = 4
          order(18,3) = 2
          order(18,4) = 1

          order(19,1) = 4
          order(19,2) = 1
          order(19,3) = 2
          order(19,4) = 3
          order(20,1) = 4
          order(20,2) = 1
          order(20,3) = 3
          order(20,4) = 2
          order(21,1) = 4
          order(21,2) = 2
          order(21,3) = 1
          order(21,4) = 3
          order(22,1) = 4
          order(22,2) = 2
          order(22,3) = 3
          order(22,4) = 1
          order(23,1) = 4
          order(23,2) = 3
          order(23,3) = 1
          order(23,4) = 2
          order(24,1) = 4
          order(24,2) = 2
          order(24,3) = 2
          order(24,4) = 1

      endif
      return
      end
      
