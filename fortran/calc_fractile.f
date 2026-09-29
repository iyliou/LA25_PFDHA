      subroutine Calc_frac ( haz1, endBR_wt, nEndBR, nZ, haz_frac )

      implicit none
      include "rupHaz.inc"
      integer MAX_MONTE, nEndBR
      parameter (MAX_MONTE=20000)
      real*8 Haz1(MAXZ,Max_TotalBr)
      real endBR_wt(Max_TotalBr), cum_wt(Max_TotalBr)
      integer nZ, iMonte, nMonte
      integer i, iseed, iBR, iz, i10, i50, i90
      real haz(MAX_MONTE), Y(MAX_MONTE)
      real haz_frac(MAXZ,20)
      real*8 meanHaz(MAXZ), sum1

c     compute cumulative weights
      cum_wt(1) = endBR_wt(1)
      DO i=2,nEndBR
        cum_wt(i) = cum_wt(i-1) + endBR_wt(i)
      enddo
      write (*,'( 2x,''test wts:'',f10.5)') cum_wt(nEndBr)

c     Reset end wt to avoid rounding errors
      cum_wt(nEndBr) = 1.0
      
      nMonte = 10000
      iseed = 1147
            
c     initialize meanHaz      
      do iz=1,nZ
        meanHaz(iz) = 0.
      enddo   
      
c     Loop over Disp values
      do iz=1,nZ
      
       do iMonte=1,nMonte

c       Sample FDM
        call GetRandom0 ( iseed, nEndBR, cum_wt, iBR, 1)

        haz(iMonte) = haz1(iz,iBR) 
        meanHaz(iz) = meanHaz(iz) + haz1(iz, iBR) 

       enddo

c      Sort hazard
       call SORT(haz,Y,nMonte)
       
       do i=1,19
        i10 = int( nMonte * 0.05 *(i) )
        haz_frac(iz,i) = haz(i10)
       enddo
       haz_frac(iz,20) = meanHaz(iz)/nMonte
       meanHaz(iz) = meanHaz(iz)/nMonte
       
      enddo
      
      return
      end

