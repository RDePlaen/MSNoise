      program tf_pws0
c
c Program to perform a zero-lag PWS in the tf-domain
c through the S-Transform.
c 
c Thanks for your interest in using this approach.
c I hope this program is useful to you. Please, note I
c can not take any warranties. I.e., use this program
c on your own risk. 
c Please, report bugs/improvements/... . Also, if you 
c need any help then do not hesitate to contact me.
c (I am preparing a new version which will be faster
c and which can handle larger dimensions.)
c
c Author: Martin Schimmel (schimmel@ictja.csic.es)
c
c Last modification: 
c   13/07/2010 : including kstnm= option
c   09/06/2013 : rm maxtr 
c   12/12/2013 : avoid pg stop in case no evla,evlo,stla,stlo available.
c   30/04/2015 : avoid traces with NaN on input.
c
c References:
c ===========
c
c This programs performs a PWS (phase weighted stack, Schimmel & Paulssen, 1997)
c in the time-frequency domain employing an S-Transform. This has been
c published in Schimmel & Gallart (2007). There you find also an explanation
c why the S-transform is analytic. To go back into the time domain,
c we use the inverse S-transform from Schimmel & Gallart (2005). This 
c inverse S-transform is also explained in the appendix of 
c Schimmel & Gallart (2007). An overview is given in Schimmel et al. (2011). 
c
c  Schimmel M., and H. Paulssen, Noise reduction and
c    detection of weak, coherent signals through phase
c    weighted stacks, GJI, 130, 497-505, 1997.
c  Schimmel M., and J. Gallart, The inverse S Transform in
c    filters with time-frequency localization , IEEE
c    Transactions on Signal Processing, 53 (11), 4417-4422,
c    doi:10.1109/TSP.2005.857065, 2005.
c  Schimmel M., and J. Gallart, Frequency-dependent phase
c    coherence for noise suppression in seismic array data,
c    J. Geophys. Res.,  112, B04303, doi:10.1029/2006JB004680, 2007.
c  Schimmel, M., Stutzmann, E., Gallart, J., Using instantaneous 
c    phase coherence for signal extraction from ambient noise data 
c    at a local to a global scale, Geophys. J. Int., 184, 494-506, 
c    doi: 10.1111/j.1365-246X.2010.04861.x, 2011.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Dimensions & parameters:
c-------------------------
c if changing dimensions remember to change them in all subroutines.
      parameter (nsmax=8192,nfmax=4096)
c     parameter (nsmax=16384,nfmax=8192)
c if changing dimensions remember to change them in all subroutines.
      character*80 filein
      character*18 par,par1
      character*8 kstnm
      real sig(nsmax),sig2(nsmax),rsig(nsmax)
      real twpi,costapl,f1,f2,t1,t2
      complex ctf(nsmax,nfmax+1)
      complex ctf2(nsmax,nfmax),ctfps(nsmax,nfmax)
      complex carray1(nsmax),carray2(nsmax),cdum,carg
      logical lrm,lkstnm

      twpi=8.0*atan(1.0)
      pi=4.0*atan(1.0)

c INPUTS:
c--------
c  Get number of arguments:
      narg=iargc()
      if (narg.lt.1) then
         write(*,*)' '
         write(*,*)'Dimensions: samples, frequencies:',nsmax,nfmax
         call usage
      endif
 
      call getarg (1,filein)
      if (narg.eq.1.and.filein(1:4).eq.'info') then
       call infooo
       stop
      endif

c DEFAULTS:
C----------
      cycle=2.
      lrm=.false.
      lkstnm=.false.
      wu=2.
      
c Get arguments (optional):
c--------------------------
      do iarg=2,narg
        call getarg(iarg,par)
        l=leng(par)
        if (par(1:3).eq.'wu=') then
          par1=par(4:l)
          read(par1,*)wu
        else if (par(1:6).eq.'kstnm=') then
          kstnm=par(7:l)
          lkstnm=.true.
        else if (par(1:4).eq.'cyc=') then
          par1=par(5:l)
          read(par1,*)cycle
        else if (par(1:4).eq.'info') then
          call infooo
          stop
        else if (par(1:2).eq.'rm') then
          lrm=.true.
          call usage                 
        endif
      enddo

c Zero arrays:
c-------------
c consumes time and is not needed for
c most compilers!
      do ns=1,nsmax
        sig(ns)=0.0
        sig2(ns)=0.0
        rsig(ns)=0.0
        do nf=1,nfmax
          ctf(ns,nf)=cmplx(0.0,0.0)
          ctfps(ns,nf)=cmplx(0.0,0.0)
          ctf2(ns,nf)=cmplx(0.0,0.0)
        enddo
      enddo

c read file with data:
c---------------------
      beg=0.0
      nskip=0
      open(1,file=filein,status="old")
      rewind(1)

      itr=0
      mxtr=0
c     do itr=1,maxtr
      do while (mxtr.eq.0)
        itr=itr+1
        read(1,'(A80)',end=80)filein
        call rsac1(filein,sig,nsamp,beg,dt,nsmax,nerr)
        if (itr.eq.1) then
          call getfhv('stla',stla,nerr1)
          call getfhv('stlo',stlo,nerr1)
          call getfhv('evla',evla,nerr1)
          call getfhv('evlo',evlo,nerr1)
        endif

c Some basic checking on the data:
c---------------------------------
        if (nerr.ne.0) then 
          write(*,*)' problem with trace number:',itr
          write(*,*)filein
          stop 'Error in reading sac file'
        endif
        if (nerr1.ne.0) then
         stla=12345
         stlo=12345
         evla=12345
         evlo=12345
        endif
        if (itr.gt.1) then
           if (nsamp.gt.nsmp) then
              write(*,*)' WARNING: use only ',nsmp,' samples on'//
     & ' trace ',itr
           endif
           if (abs(dt-dt1).gt.dt1*0.01) then
c            stop ' Traces have different dt ! '
             write(*,*)' ATTENTION: trace ',itr,' has a different dt !'
             write(*,*)' ATTENTION: skipping trace ',itr,' .'
             nskip=nskip+1
             goto 999
           endif
           if (abs(beg1-beg).gt.dt1) then
c            stop ' Traces have different begin time.'
             write(*,*)' ATTENTION: trace ',itr,' has a different beg !'
             write(*,*)' ATTENTION: skipping trace ',itr,' .'
             nskip=nskip+1
             goto 999
           endif
        endif
        call getfhv('depmax',depmax,nerr1)
        if (isnan(depmax)) then
          write(*,*)' ATTENTION: trace ',itr,' with problems (NaNs)'
          write(*,*)' ATTENTION: skipping trace ',itr,' .'
          nskip=nskip+1
          goto 999
        endif

C RM MEAN
c--------
        if (lrm) then
          rmean=0.0
          do i=1,nsamp
            rmean=rmean+sig(i)
          enddo
          do i=1,nsamp
            sig(i)=sig(i)-rmean
          enddo
        endif
c linear t-domain stack
c-----------------------
        do i=1,nsamp
           rsig(i)=sig(i)+rsig(i)
        enddo

c get power for FFT:
c--------------------
        if (itr.eq.1) then
          beg1=beg
          dt1=dt
          npow=1
          nsmp=2
 5        if (nsamp.le.nsmp) goto 7
          npow=npow+1
          nsmp=nsmp*2
          goto 5
          stop ' problem with power?'
 7        nsmp2=nsmp/2
c         df=1./(dt*float(nsmp))
        endif
        
c S-transform
c------------
        call st_4ward(sig,nsamp,cycle,dt,ctf,nsmp,npow)

c tf-domain stack ctf2 & PS
c-----------------------
        do ns=1,nsmp
          do nf=1,nsmp2
             cdum=ctf(ns,nf)
             ctf2(ns,nf)=ctf2(ns,nf)+cdum
             rdum=cabs(cdum)
             if (rdum.lt.0.00001) then
               cdum=cmplx(0.,0.)
             else
               cdum=cdum/(cabs(cdum))
             endif
             ctfps(ns,nf)=ctfps(ns,nf)+cdum
          enddo
        enddo
        
 999    continue
        if (itr.le.3) write(*,*)' finished with trace #:',itr
        if (itr.eq.4) write(*,*)'           ...'
        if (mod(itr,20).eq.0) then
           write(*,*)' finished with trace #:',itr
           write(*,*)'           ...'
        endif
c end loop while over all traces!
      enddo
 80   continue
      mtr=itr-1-nskip
      write(*,*)' total # traces       :',mtr

c tf-PWS:
c--------
      rmtr=float(mtr)
      do ns=1,nsmp
        do nf=1,nsmp2
           cdum=ctf2(ns,nf)/rmtr
           rdum=cabs(ctfps(ns,nf))/rmtr
           ctf(ns,nf)=cdum*(rdum**wu)
        enddo
      enddo

c INVERSE TRANSFORM:
c-------------------
      write(*,*)'   ... just before inv transform.'
      cdum=cmplx(0.0,-twpi/float(nsmp))
c     rdum=float(nsmp)*twpi*2.
      rdum=float(nsmp)*sqrt(twpi)
      do ntime=1,nsamp
         carg=cdum*(ntime-1)
         do nfre=1,nsmp2
            carray1(nfre)=ctf(ntime,nfre)*rdum
         enddo
         sig(ntime)=diftc_fast_f(carray1,nsmp,nsmp2,carg)
      enddo

c Normalize LS:
      do ntime=1,nsamp
         rsig(ntime)=rsig(ntime)/rmtr 
      enddo

      write(*,*)'   ... write sac traces: tl.sac, tf-pws.sac'
      if (lkstnm) then
      call wrsac('tl.sac',kstnm,beg,dt,rsig,nsamp,evla,evlo,stla,stlo)
      call wrsac('tf-pws.sac',kstnm,beg,dt,sig,nsamp,evla,evlo,
     & stla,stlo)
      else
      call wrsac('tl.sac','t-lin',beg,dt,rsig,nsamp,evla,evlo,stla,stlo)
      call wrsac('tf-pws.sac','tf-pws',beg,dt,sig,nsamp,evla,evlo,
     & stla,stlo)
      endif

      end
c--------------------------------------------------------------
      SUBROUTINE usage
      write(*,*)'      '
      write(*,*)'USAGE:  tf_pws0 filename parameter_list '
      write(*,*)'      '
      write(*,*)'filename: One file name (SAC) per line; '
      write(*,*)'          traces must have same b,nsmpl&dt.'
      write(*,*)'      '
      write(*,*)'parameter_list: rm cyc= wu= info kstnm='
      write(*,*)'                 '
      write(*,*)'parameters are optional and can be provided in'
      write(*,*)'arbitrary order without any blank around =.'
      write(*,*)' '
      write(*,*)'rm     : remove mean from input data. '
      write(*,*)'cyc=   : num. of cycles in 2*std (default cyc=2)'
      write(*,*)'wu=    : PWS power (default: wu=2)'
      write(*,*)'kstnm= : optionally specify sac-header (char*8)'
      write(*,*)'info   : write background and main references'
      write(*,*)'         to screen. Just type: tf_pws0 info'
      write(*,*)
      write(*,*)'DEFAULTS: cyc=2 wu=2'
      write(*,*)' '
      write(*,*)'AUTHOR: Martin Schimmel, 13/07/2010'
      write(*,*)' '
      write(*,*)'Please, do not hesitate to send bugs, comments'//
     & ', improvements,'
      write(*,*)'nice results, ... to schimmel@ictja.csic.es .'
      write(*,*)' '
      stop
      return
      end
c--------------------------------------------------------------
      SUBROUTINE st_4ward(sig,nsamp,cycle,dt,ctf,nsmp,npow)
c
c INPUT: sig,nsamp,cycle,dt,nsmp,npow
c * sig(nsamp) is real time series with sample interval dt.
c * cycle is the number of periods that equals the 2 std of
c   the Gaussian window.
c * nsmp=2**npow  (==> for FFT)
c * dt is sample interval (time).
c OUTPUT: ctf
c * ctf(nsamp,nsmp2+1) is the complex t-f representation of
c   sig(nsamp).
c   nsamp in ctf is the center sample of 'moving Gauss-window'.
c
c  schimmel@ija.csic.es,  07.05.2004
c  last change 22.08.06
c=================================================================

c     parameter (nsmax=16384,nsctf=16384,nfctf=8193)
      parameter (nsmax=8192,nsctf=8192,nfctf=4097)
c     parameter (nsmax=4096,nfmax=2048)
c nfmax=nsmax/2 + 1, nsctf=2**n, nsctf>nsmax
      real sig(nsmax)
      complex csig(nsctf),ctf(nsctf,nfctf)
      complex cdum(nsctf)
      integer npow,nsmp,nsamp
      real pi,fact,cycle

      pi=4.0*atan(1.0)
      fact=cycle/2.
      nsmp2=nsmp/2

      do ns=1,nsamp
        csig(ns)=cmplx(sig(ns),0.0)
      enddo
c to go the safe way:
      do ns=nsamp+1,nsctf
        csig(ns)=cmplx(0.0,0.0)
      enddo
c
c Perform FFT:
c-------------
      rn=1.
      call clogc(npow,csig,rn,dt)

c=======================
c Start frequency loop:
c=======================
      pp=-2.0*pi*pi*fact*fact

c for n=0 (zero frequency nf=1):
c (zero mean t-series assumed!!!!)
c---------------------------------
      do ns=1,nsamp
        ctf(ns,1)=cmplx(0.0,0.0)
      enddo

c and now remaining positive frequencies:
c----------------------------------------
      do nf=2,nsmp2+1
        nf1=nf-1
c zero negative freqs:
        do mf=1,nsmp
          cdum(mf)=cmplx(0.0,0.0)
        enddo
        rpp=pp/float(nf1)/float(nf1)
c       do mf=1,nsmp+1
        do mf=1,nsmp
          rmf=float(mf-1)
          rns=float(nsmp-mf+1)
          argp=rpp*rmf*rmf
          argn=rpp*rns*rns
          if (argp.gt.-25.or.argn.gt.-25) then
c shift Fourier spectrum csig(mf) by freq nf-1 to csig(mf+nf-1):
            mn=mf+nf1
            if (mn.gt.nsmp) mn=mn-nsmp
c multiply with localizing Gaussian (generalized voice Gaussian):
            rr=exp(argp)+exp(argn)
            cdum(mf)=csig(mn)*rr
          endif
        enddo
         
c perform ifft to obtain temporal spectral localizations:
        rn=-1.
        call clogc(npow,cdum,rn,dt)
c store in matrix:
        do ns=1,nsamp
          ctf(ns,nf)=cdum(ns)
        enddo
      enddo

      return
      end
c----------------------------------------------------------------
      SUBROUTINE infooo

      write(*,*)' '
      write(*,*)'tf_pws0.f: Zero-lag PWS in the tf-domain '
      write(*,*)' '
      write(*,*)'This program performs a PWS (phase weighted '//
     & 'stack, Schimmel & '
      write(*,*)'Paulssen, 1997) in the time-frequency domain'//
     & ' employing an'
      write(*,*)'S-Transform. This has been published in'//
     & ' Schimmel & Gallart (2007).'
      write(*,*)'There you find also when and why the S-transform'//
     & ' is analytic. To go '
      write(*,*)'back into the time domain, we use the inverse'//
     & ' S-transform from' 
      write(*,*)'Schimmel & Gallart (2005). This inv. S-transform '//
     & 'is also explained '
      write(*,*)'in the appendix of Schimmel & Gallart (2007). '//
     & 'Schimmel et al. (2011)'
      write(*,*)'show how these straegies can be used '//
     & 'with ambient seismic noise.'
      write(*,*)
      write(*,*)'Informations on how to discretise the S-transform'//
     & ' without side effects'
      write(*,*)'are published in Simon et al. (2007). Further,'//
     & ' Ventosa et al. (2008)'
      write(*,*)'show the relations between the S-transform and'//
     & ' wavelet transform.'
      write(*,*)
      write(*,*)'   Schimmel M., and H. Paulssen, Noise reduction '//
     & 'and'
      write(*,*)'     detection of weak, coherent signals through '
      write(*,*)'     phase weighted stacks, GJI, 130, 497-505, 1997.'
      write(*,*)'   Schimmel M., and J. Gallart, The inverse S '//
     & 'Transform in'
      write(*,*)'     filters with time-frequency localization , '//
     & 'IEEE'
      write(*,*)'     Transactions on Signal Processing, 53 (11), '//
     & '4417-4422,'
      write(*,*)'     doi:10.1109/TSP.2005.857065, 2005.'
      write(*,*)'   Schimmel M., and J. Gallart, Frequency-dependent '//
     & 'phase'
      write(*,*)'     coherence for noise suppression in seismic '//
     & 'array data,'
      write(*,*)'     J. Geophys. Res., 112, B04303, '//
     & 'doi:10.1029/2006JB004680, 2007.'
      write(*,*)'   Schimmel, M., Stutzmann, E., Gallart, J., Using '//
     & 'instantaneous '
      write(*,*)'     phase coherence for signal extraction from '//
     & 'ambient noise'
      write(*,*)'     data at a local to a global scale, Geophys. J. '//
     & 'Int., 184, '
      write(*,*)'     494-506, doi:10.1111/j.1365-246X.2010.04861.x, '//
     & '2011.'
      write(*,*)'   Simon et al., The S-transform and its inverses:'//
     & ' side effects '
      write(*,*)'     of discretising and filtering, IEEE '//
     & 'Transactions on Signal'
      write(*,*)'     Processing, 55, 4928-4937, '//
     & 'doi:10.1109/TSP.2007.897893, 2007.'
      write(*,*)'   Ventosa et al., S-transform from a wavelets point'//
     & ' of view,'
      write(*,*)'     IEEE Transactions on Signal Processing, 56, '//
     & '2771-2780,'
      write(*,*)'     doi:10.1109/TSP.2008.917029, 2008.'
      write(*,*)
      write(*,*)' Author: Martin Schimmel (schimmel@ija.csic.es)'
      write(*,*)
      write(*,*)' Last modification: 09/06/2013'
      write(*,*)

      return
      end
c--------------------------------------------------------------------
      SUBROUTINE wrsac(name,kst,b,dt,y,nsmpl,evla,evlo,stla,stlo)
 
      real y(1),b,dt,gcarc,evdp,baz
      integer nsmpl
      character*(*) kst,name
 
      call newhdr
      call setnhv('npts',nsmpl,nerr)
      call setfhv('delta',dt,nerr)
ccc   call setfhv('gcarc',gcarc,nerr)
ccc   call setfhv('evdp',evdp,nerr)
ccc   call setfhv('baz',baz,nerr)
      call setkhv('kstnm',kst,nerr)
c     call setfhv('t7',t1,nerr)
      call setfhv('evla',evla,nerr)
      call setfhv('evlo',evlo,nerr)
      call setfhv('stla',stla,nerr)
      call setfhv('stlo',stlo,nerr)
 
      e=b+(nsmpl-1)*dt
      o=0.
 
c     call setfhv('e',e,nerr)
c     call setfhv('o',o,nerr)
      call setfhv('b',b,nerr)
c     call setihv('iztype','io',nerr)
c
      call wsac0(name,dum,y,nerr)
      if(nerr.ne.0) stop 'Error in writing output file'
 
      return
      end
c--------------------------------------------------------------
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
c--- the second half of x contains the spectrum for negative frequencies
c--- such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
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
c-------------------------------------------------------------
      integer FUNCTION leng(char)
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
c-------------------------------------------------------------
        REAL FUNCTION diftc_fast_f(y,nsmp,nsmp2,carg)
c
c Discrete Inverse Fourier Transformation / |f|:
c (Note: constant factors are not included.)
c
c  y        is array with complex spectrum (positive freqs).
c  nsmp     is # of samples in time domain.
c  nsmp2    is the length of data vector y (pos. freqs).
c  ntime    is desired time index.
c  carg     is cmplx(0.0,twpi*ntime/nsmp)
c===============================================================
        complex y(1),carg,yy,cdum,cphs1,cphs2
        integer nsmp2,nsmp

        cdum=cmplx(0.D0,0.D0)
        do nm=1,nsmp2
           yy=y(nm)/float(nm)
           cphs1=yy*cexp(carg*(nm-1))
           cphs2=conjg(yy)*cexp(carg*(nsmp-nm+1))
           cdum=cdum+cphs1+cphs2
        enddo
        diftc_fast_f=real(cdum)/float(nsmp)

       end
c---------------------------------------------------------------
