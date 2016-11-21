      program norm_tavg_white
c
c  pg for 1-bit normalization and spectral whitening. 
c  Use on own risk!
c
c  (schimmel@ictja.csic.es)
c----------------------------------------------------------

c     parameter (max=9999999,nsmax=16384)
      parameter (max=9999999,nsmax=1048576)
      character*80 name
      character*15 par,par1
      real z(max), zz(max)
      logical lavg,lwit
      complex csig(nsmax),carg,ccarg

cccccccccccccccccccccc
c get the input files:
cccccccccccccccccccccc
      lwit=.false.
      lavg=.false.

c number of arguments: narg
      narg=iargc()
      rrr=-1
 
      if(narg.lt.1) call usage
 
      call getarg(1,name)
      do iarg=2,narg
        call getarg(iarg,par)
        l=leng(par)
        if (par(1:2).eq.'N=') then
          par1=par(3:l)
          read(par1,*)nnn
          lavg=.true.
        elseif (par(1:2).eq.'T=') then
          par1=par(3:l)
          read(par1,*)rrr
          lavg=.true.
        elseif (par(1:5).eq.'white') then
          lwit=.true.
        else
          call usage
        endif
      enddo

ccccccccccccccccc
c read trace:
ccccccccccccccccc
       call rsac1(name,z,nlen,beg1,dt1,max,nerr)
       if (nlen.gt.max) stop ' change dimensions'

       if (rrr.gt.0)  nnn=nint((rrr/2*dt1))
cccccccccccccc
c time window:
cccccccccccccc
cc    n1=nint((t1-beg1)/dt1) + 1
cc    n2=nint((t2-beg1)/dt1) + 1
cc    if (n2.gt.nlen) stop ' Is T2 too large? '

ccccccccccc
c     norm:
ccccccccccc
      if (lavg) then
      do i=1,nlen
        m1=i-nnn
        m2=i+nnn
        if (m1.lt.1) m1=1
        if (m2.gt.nlen) m2=nlen
        mmm=m2-m1+1
        rmean=0.
        do m=m1,m2
          rmean=rmean+abs(z(m))
        enddo
        if (rmean.gt.0.) then
          zz(i)=z(i)*float(mmm)/rmean
        else
          zz(i)=z(i)
        endif
      enddo
      endif

ccccccccccc
c   whiten:
ccccccccccc
      if (lwit) then
       if (nlen.gt.nsmax) stop 'change nsmax for whitening.'
       if (lavg) then
        do i=1,nlen
         csig(i)=cmplx(zz(i),0.)
        enddo
       else
        do i=1,nlen
         csig(i)=cmplx(z(i),0.)
        enddo
       endif
       do i=nlen+1,nsmax
        csig(i)=cmplx(0.,0.)
       enddo
c get power
       npow=2
  5    npow=npow+1
       if (2**npow.lt.nlen) goto 5
       nsmp=2**npow
       nsmp2=nsmp/2
c FFT
       rn=1.0
       dt=1.
       call clogc(npow,csig,rn,dt)
c whiten
       do n=1,nsmp2+1
         rr=cabs(csig(n))
         if (rr.gt.0.0001) then
           csig(n)=csig(n)/rr
         else
           rr=atan2(aimag(csig(n)),real(csig(n)))      
           csig(n)=cexp(cmplx(0.,rr))
         endif
       enddo
       do n=nsmp2+2,nsmp
        nn=nsmp-n+2
        csig(n)=conjg(csig(nn))
       enddo
      csig(nsmp2+1)=0.
c IFFT:
       rn=-1.
       call clogc(npow,csig,rn,dt)

       do i=1,nlen
        zz(i)=real(csig(i))
       enddo

      endif

ccccccccccccccccccc
c write sac file:
ccccccccccccccccccc
      l=leng(name)
      name(1:l+1)=name(1:l)//"n"
      call wsac0(name(1:l+1),xdum,zz,nerr) 

      end
c______________________________________________________________________
      integer function leng(char)
c     leng is the position of the last non blank character
c     in the character variable char
      character char*(*)
      integer*4 l
      l=len(char)
      do 1 leng=l,1,-1
    1 if (char(leng:leng) .ne. ' ') return
      leng=0
      return
      end
c______________________________________________________________________
      SUBROUTINE usage

      write(*,*)
      write(*,*)'running-absolute-mean-normalization and/or spectral '//
     & 'whitening'
      write(*,*)
      write(*,*)'USAGE: norm_tavg_white inputfile paramters'
      write(*,*)' ' 
      write(*,*)'inputfile: sac_file '
      write(*,*)'   '
      write(*,*)'  N=     : 2N+1 is the window size of the mov window.'
      write(*,*)'           s_new_j = s_j / w_j with '
      write(*,*)'           w_j = ( sum_k |s_k| ) / ( 2N + 1),  k'//
     & ' in [j-N,j+N]'
      write(*,*)'           N about half maximum period (Bensen et'//
     & ' al. 2007).'
      write(*,*)'           N=0 ==> 1-bit normalization'
      write(*,*)'  T=     : (alternative to N=) N=nint((T/2dt))'
      write(*,*)'  white  : Spectral whitening through FFT.'
      write(*,*)'           Whitening is applied after running norm.'
      write(*,*)
      write(*,*)'output file: program appends n to input file name.'
      write(*,*)
      write(*,*)' EXAMPLE 1 (1bit normalization):'
      write(*,*)'   norm_tavg_white my_tseries.sac N=0 '
      write(*,*)'   output: my_tseries.sacn '
      write(*,*)' EXAMPLE 2 (spectral whitening):'
      write(*,*)'   norm_tavg_white my_tseries.sac white '
      write(*,*)'   output: my_tseries.sacn '
      write(*,*)' EXAMPLE 3 (1bit + spectral whitening):'
      write(*,*)'   norm_tavg_white my_tseries.sac N=0 white '
      write(*,*)'   output: my_tseries.sacn '
      write(*,*)' '
      write(*,*)' Bugs, questions, comments: schimmel@ictja.csic.es'
      write(*,*)' '
      write(*,*)' LAST CHANGE: 11.07.10'
      stop
      return
      end
c-------------------------------------------------------------------
      SUBROUTINE CLOGC(n,x,sigm,dt)

c--- performs fft on signals with length 2**n and sampling interval
c--- of dt seconds (if in the time domain; notice that dt*df=1/2**n).
c--- the signal is stored in x. it may be complex.
c--- the spectrum is returned in x. it is almost always complex.
c--- a time-to-frequency transform is done with sign=+1. (conform
c--- the convention adopted in aki and richards - the alternative
c--- convention may be obtained by taking complex conjugates after
c--- the call to clogc).
c--- the normalization factor 1./twopi occurs in the frequency-to
c--- time transform (again aki&richards).
c--- normalization is such that physical dimensions are respected.
c--- thus, if the time signal is dimensioned in meters, the
c--- resulting spectral density in x is in meters/hz. for example,
c--- if the time signal is the unit sinc function of width dt, centered
c--- at t=0, the spectral density is dt for all values of the frequency.c
c--- array locations: if x contains the spectrum, it has the spectrum
c--- for positive frequencies in the first 2**n/2+1 elements, such that
c--- x(1) is at 0 hz, x(2) at df hertz, and x(2**n/2+1) at the nyquist,
c--- where df=1./(2**n*dt) and the nyquist is 1./(2*dt) hz.
c--- the second half of x contains the spectrum for negative frequenciesc--- such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
c--- if x contains the time signal, x(1) is at time 0, x(2)
c--- at time dt etc.
c
      dimension x(1),m(25)
      complex x,wk,hold,q
      if(sigm.ge.0.) then
        sign=1.
      else
        sign=-1.
      endif
      lx=2**n
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=sign*6.283185308*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q
    2 continue
      do 3 i=2,n
      ii=i
      if(k.lt.m(i)) go to 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if(k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if(k.lt.m(i)) go to 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if(sign.gt.0.) go to 9
      flx=flx*dt
      do 8 i=1,lx
    8 x(i)=x(i)/flx
      return
    9 do 10 i=1,lx
   10 x(i)=x(i)*dt
      return
      end
c--------------------------------------------------------------------
