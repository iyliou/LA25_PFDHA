
      subroutine GetRandom0 ( iseed, n, wt, iSave, iflag )

      implicit none

      integer iseed, i1, isave, n1, n, i, iflag, iflt
      real ran1, wt(1), x
      
c     Get random number
      x = ran1( iseed )

      do i=1,n
        if ( x .le. wt(i) ) then
          iSave = i
          return
        endif
      enddo
      
      write (*,*) ' Get Random Number 0'
      write (*,'( 2x,''Error - bad ran number or weights'')')
      write (*,*) ' Random Number        = ', x
      write (*,*) ' Fixed Parameter      = ', i1
      write (*,*) ' Number of Parameters = ', n
      write (*,*) 'iflag=',iflag
      do i=1,n
         write (*,'(2e12.3)') wt(i), x
      enddo
      stop 99
      end
      
c ----------------------------

c ----------------------------

      function Ran1 ( idum )

      implicit none

c     Random number generator, From numerical recipes
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807, im=2147483647, am=1./im,iq=127773,ir=2836,
     1      ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      
      if (idum .le. 0 .or. iy .eq. 0 ) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if( idum .lt. 0) idum=idum+im
          if (j .le. ntab) iv(j)=idum
        enddo
        iy = iv(1)
      endif
      k = idum/iq
      idum=ia*(idum-k*iq) - ir*k
      if ( idum .lt. 0 ) idum = idum + im
      j = 1 + iy/ndiv
      iy = iv(j)
      iv (j) = idum
      ran1 = min(am*iy,rnmx)
      return
      end

