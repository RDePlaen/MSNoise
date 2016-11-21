      program pcc5red
c
c  Phase Cross-Correlation! 
c
c  Program is distributed for research purposes and in the hope 
c  that it will be useful however without any warranties. 
c  Use it on your own risk. Please, report bugs/improvements/... .  
c
c  This program performs a phase cross-correlatation (PCC)
c  of two t-series in the time domain.
c  The 1st time series contains the 'pilot wavelet'. 
c  PCC is computed also for the partial overlaps. 
c  The begin and end lag (units: samples or seconds) can be specified.
c  A positive lag means that the 1st trace is shifted
c  by a positive lag time with respect to the 2nd trace.
c
c  This pg requires input traces with same starting time!!!!
c  Cut the traces with sac or consult pcc3 or pcc4 if you 
c  have sac traces with different start time. Up to now 
c  only pcc5 handels asc & bin i/o.
c
c  This is a derivative of pcc4, several changes have been done.
c  Simplifications are required for a straight final code. Many 
c  issues are solved in a complicated way. Essentially, since this 
c  pg has been derived from the derivative of the derivative of an
c  ancient test-pg based on a different philosophy. 
c
c  The pg should be coded newly for more efficiency.
c
c  Last change:   08.07.02 merge ccgn_t & pcc
c  Last change:   05.11.10 speeding up, including nll  
c                          (env has not yet tested with nll).
c  Last change:   12.12.11 speeding up following ideas provided by
c                          Sergi Ventosa.
c  Last change:   14.12.11 including pow in PCC norm.
c  Last change:   21.03.13 pcc5d ==> pcc5e
c                          Permit to compute ccgn alone.
c                          Rm need to provide nlen1, nlen2
c  Last change:   29.05.13 pcc5e ==> pcc5f
c                          Analytic signal normalization: A possible 
c                          division by zero may occurr (should not 
c                          happen for real data) and cause an NaN output.
c  Last change:   05.03.14 pcc5f ==> pcc5g
c                          replace ccgn routine.
c  NEXT changes:  interpolation, different beg time, faster fft,...
c  Reference:     Schimmel, BSSA, 1999.     
c  Author:        Martin Schimmel (schimmel@ictja.csic.es)
c
c  nlen1,nlen2  # of samples of 1st and 2nd t-series.
c  nl1,nl2      min, max lag (samples)
c  n11,n12      1st and last samples of pilot on t-series.
c  n21,n22      1st and last samples of trace on t-series.
c--------------------------------------------------------------
c     parameter (max=8192)
c     parameter (max=32768)
c     parameter (max=65536)
c     parameter (max=131072)
c     parameter (max=262144)
c     parameter (max=524288)
c     parameter (max=1048576)
      parameter (max=8388608)
      complex csig1(max),csig2(max),cdum
      real sig1(max),sig2(max),rcc(max),wu
      complex pilot(max),trace(max)
      character*80 name1, name2 
      character*14 nameout1t, nameout1, nameout2
      character kstnm*8, nameout3*12
c I/O format:
      logical lisac,liasc,libin
      logical losac,loasc,lobin
      logical lccgn,lenv,lnn,lpcc
      logical lpow
c
      twopi=6.28318530717958
      pi=twopi/2.
      mmax=max

cccccccccccccc
c Zero arrays:
cccccccccccccc
      do ns=1,max
         sig1(ns)=0.
         sig2(ns)=0.
      enddo

ccccccccccccccccccccccc
c get input parameters:
ccccccccccccccccccccccc
      call getpar(name1,name2,dt,wu,liasc,lisac,libin,
     &  loasc,losac,lobin,lpcc,lccgn,lenv,lnn,lpow,
     &  nlen1,nlen2,n11,n12,t11,t12,nl1,nl2,tl1,tl2,nll,mmax)

ccccccccccccccccc
c get input data:
ccccccccccccccccc
      if (lisac) then
	write(*,*)'THIS IS A SAC FILE!!'
        call sactraces(name1,name2,sig1,sig2,nlen1,nlen2,
     &     max,dt2,beg,gcarc,evdp,kstnm,baz,s1la,s1lo,s2la,s2lo,lnn)
        if (dt.ne.dt2) then
           dt=dt2
           write(*,*)'use dt from sac traces!! dt= :',dt
        endif
      write(*,*)'length of t-series (samples) : ',nlen1,nlen2
      else if (libin) then
        open(unit=1,form='UNFORMATTED',file=name1,err=825)
        open(unit=2,form='UNFORMATTED',file=name2,err=826)
        goto 830
825     write(nstderr,'(a)') 'Can''t open file : ',name1
        stop
826     write(nstderr,'(a)') 'Can''t open file : ',name2
        stop
830     rewind 1
        rewind 2
        read(1,err=999)(sig1(j),j=1,nlen1)
        read(2,err=999)(sig2(j),j=1,nlen2)
        close(unit=1)
        close(unit=2)
      else if (liasc) then
        open(unit=1,file=name1,err=835)
        open(unit=2,file=name2,err=836)
        goto 840
835     write(nstderr,'(a)') 'Can''t open file : ',name1 
        stop
836     write(nstderr,'(a)') 'Can''t open file : ',name2 
        stop
840     rewind 1
        rewind 2
        read(1,*,err=999)(sig1(j),j=1,nlen1)
        read(2,*,err=999)(sig2(j),j=1,nlen2)
        close(unit=1)
        close(unit=2)
      else 
         call usepcc(mmax)
      endif

      if (n11.eq.0.and.t11.eq.0.) n11=1
      if (n12.eq.0.and.t12.eq.0.) n12=nlen1
      if (n21.eq.0) n21=1
      if (n22.eq.0) n22=nlen2
      if (nll.lt.1) nll=1
      if (tl1.ne.0) nl1=nint(tl1/dt)
      if (tl2.ne.0) nl2=nint(tl2/dt)
      if (t11.ne.0) n11=nint((t11-beg)/dt)+1
      if (t12.ne.0) n12=nint((t12-beg)/dt)+1
      write(*,*)'pilot n11,n12= (samples)     :',n11,n12
      write(*,*)'sample lag nl1,nl2=          :',nl1,nl2
      write(*,*)

ccccccccccccccccccccccccccccccc
c CCGN TIME DOMAIN computation:
ccccccccccccccccccccccccccccccc
      if (lccgn) then
         call ccgn_t(sig1,sig2,rcc,nlen2,nccgn,
     &      n11,n12,nl1,nl2,nll)
      endif

cccccccccccccccccccccccccccccc
c PCC TIME DOMAIN computation:
cccccccccccccccccccccccccccccc
      if (lpcc) then
c zero arrays:
        do ns=1,max
          pilot(ns)=cmplx(0.,0.)
          trace(ns)=cmplx(0.,0.)
          csig1(ns)=cmplx(0.,0.)
          csig2(ns)=cmplx(0.,0.)
       enddo
ccccccccccccccccccccccc
c get ANALYTIC signals:
ccccccccccccccccccccccc
        call analytic(sig1,csig1,nlen1)
        call analytic(sig2,csig2,nlen2)

cccccccccccccccccccccc
c get INSTANT. PHASES: 
cccccccccccccccccccccc
c & Xtract pilot:
        n1dum=0
        do ns=n11,n12
           n=ns-n11+1
           cdum=csig1(ns)
           rdum=cabs(cdum)
           if (rdum.lt.0.0000000001) then
             pilot(n)=cmplx(0.,1.)
             n1dum=n1dum+1
           else
             pilot(n)=cdum/rdum
           endif
        enddo
c Xtract trace:
c n21-n22 is not necessary since specified through nl1 & nl2!
        n2dum=0
        do ns=1,nlen2
           cdum=csig2(ns)
           rdum=cabs(cdum)
           if (rdum.lt.0.0000000001) then
             n2dum=n2dum+1
             trace(ns)=cmplx(1.,0.)
           else
             trace(ns)=cdum/rdum
           endif
        enddo
        npilot=n12-n11+1
        ntrace=nlen2

        if (n1dum.ne.0) then
         write(*,*)' Warning: analytic trace (1) with 0 '//
     &   'amplitude samples:',n1dum
        endif
        if (n2dum.ne.0) then
         write(*,*)' Warning: analytic trace (2) with 0 '//
     &   'amplitude samples:',n2dum
        endif

c     ll=0
        nlag1=n11+nl1
        nlag2=n11+nl2
        nnn=nlag2-nlag1+1
        if (nnn.gt.max) then
          nlag2=nlag1+max-1
          write(*,*)' NOT ALL LAGS POSSIBLE !'
        endif
           
        if (lpow) then
         call pcc_time_exp_out3(sig1,pilot,trace,nlag1,nlag2,nll,
     &                        wu,npilot,ntrace,ll,n11)
        else
         call pcc_time_exp_in3(sig1,pilot,trace,nlag1,nlag2,nll,
     &                        wu,npilot,ntrace,ll,n11)
        endif

        nsmpl=ll

c PCC output:
c=============
        dtt=nll*dt
        if (lenv) then
          call envel(sig1,sig2,nsmpl)
          if (loasc) then
            tb=nl1*dt      
            open(1,file='pcc_env.asc')
            do ns=1,nsmpl
              write(1,*)tb+(ns-1)*dtt,sig2(ns)
            enddo
            close(1)
          endif
          if (lobin) then
            open(1,file='pcc_env.bin',form='UNFORMATTED')   
            write(1)(sig2(ns),ns=1,nsmpl)
            close(1)
          endif
          if (losac) then
            tb=nl1*dt      
            kstnm="pcc"
            call wrsac('pcc_env.sac',kstnm,tb,dtt,gcarc,
     &        evdp,sig2,nsmpl,lnn,s1la,s1lo,s2la,s2lo)
          endif
        endif

        if (loasc) then
          tb=nl1*dt      
          open(1,file='pcc.asc')
          do ns=1,nsmpl
            write(1,*)tb+(ns-1)*dtt,sig1(ns)
          enddo
          close(1)
        endif
        if (lobin) then
          open(1,file='pcc.bin',form='UNFORMATTED')   
          write(1)(sig1(ns),ns=1,nsmpl)
          close(1)
        endif
        if (losac) then
          tb=nl1*dt      
          kstnm="pcc"
          call wrsac('pcc.sac',kstnm,tb,dtt,gcarc,evdp,sig1,nsmpl,
     &    lnn,s1la,s1lo,s2la,s2lo)
        endif

      endif

c CCGN output:
c=============
      if (lccgn) then
        dtt=nll*dt
        if (lenv) then
          call envel(rcc,sig2,nccgn)
          if (loasc) then
            tb=nl1*dt      
            open(1,file='ccgn_env.asc')
            do ns=1,nccgn
              write(1,*)tb+(ns-1)*dtt,sig2(ns)
            enddo
            close(1)
          endif
          if (lobin) then
            open(1,file='ccgn_env.bin',form='UNFORMATTED')   
            write(1)(sig2(ns),ns=1,nccgn)
            close(1)
          endif
          if (losac) then
            tb=nl1*dt      
            kstnm="ccgn"
            call wrsac('ccgn_env.sac',kstnm,tb,dtt,gcarc,
     &          evdp,sig2,nccgn,lnn,s1la,s1lo,s2la,s2lo)
          endif
        endif
        if (loasc) then
          tb=nl1*dt      
          open(1,file='ccgn.asc')
          do ns=1,nccgn
            nns=ns-1
            write(1,*)tb+float(nns)*dtt,rcc(ns)
          enddo
          close(1)
        endif
        if (lobin) then
          open(1,file='ccggn.bin',form='UNFORMATTED')   
          write(1)(rcc(ns),ns=1,nccgn)
          close(1)
        endif
        if (losac) then
          tb=nl1*dt      
          kstnm="ccgn"
          call wrsac('ccgn.sac',kstnm,tb,dtt,gcarc,
     &       evdp,rcc,nccgn,lnn,s1la,s1lo,s2la,s2lo)
        endif
      endif

      goto 990

 999  stop 'Problem with bin data.'
 990  continue
      end
c==============================================================
c                 S U B R O U T I N E S
c==============================================================
       SUBROUTINE ccgn_t(x,y,z,ny,nz,nx1,nx2,nl1,nl2,nll)
c
c Time-domain cross-correlation with geometrical normalization!
c (Partial overlap is permitted!)
c
c  INPUT:
c  x(t),y(t), pilot & t-series with nsamp samples.
c  The pilot is isolated from the t-series x by nx1-nx2.
c  nl1,nl2 define the relative lag times. 
c  (All units are sample index numbers) 
c
c  OUTPUT:
c  z(t)=Cxy(t)/{sqrt[Cxx(t)*Cyy(t)]}
c
c  Cxy(t)=int{x(tau)y(tau+t)d_tau}
c  Cxx(t)=int{x(tau)**2 d_tau}
c  Cyy(t)=int{y(tau+t)**2 d_tau}
c
c  References:  e.g. Clayton et al., 1976 BSSA 66, p 325-326,
c               Schimmel, BSSA, 1999.  
c 2013 Changed!!!
c--------------------------------------------------------------
c     parameter (maxns=8192)
c     parameter (maxns=32768)
c     parameter (maxns=65536)
c     parameter (maxns=131072)
c     parameter (maxns=262144)
c     parameter (maxns=524288)
c     parameter (maxns=1048576)
      parameter (maxns=8388608)
      real x(maxns),y(maxns),z(maxns)
      real cxxyy(maxns),ccsig(maxns)
      integer ny,nz,nx1,nx2,nl1,nl2,npilot

      npilot=nx2-nx1+1
      if (nx2.gt.maxns) stop 'Problem with pilot dimensions.'
      if (npilot.gt.maxns) stop 'Change pilot or pg dimensions.'
      if (ny.gt.maxns) stop 'Change trace or pg dimensions.'

c  Set proper rlim! This is only required for noise free data.
c  rlim protects the normalization with very small numbers.
c  Smallest value: rlim * rnorm_max
      rlim=.00001

cccccccccccccc
c ZERO arrays:
cccccccccccccc
      do ns=1,maxns
        z(ns)=0.
        cxxyy(ns)=0.
        ccsig(ns)=0.
      enddo

ccccccccccccccccc
c Get cross-corr:
ccccccccccccccccc
      n1=nx1+nl1
      n2=nx1+nl2
      if (n1.lt.1-npilot) then
        n1=2-npilot
        nl1=n1-nx1
        write(*,*)' Minimum lag has been changed: Check it!'
        write(*,*)' new nl1=',nl1
      endif
      if (n2.gt.ny) then 
        n2=ny
        nl2=n2-nx1
        write(*,*)' Maximum lag has been changed: Check it!'
        write(*,*)' new nl2=',nl2
      endif

c     nz=n2-n1+1
c i is the sample index where to position the 1st sample
c of the pilot trace.
      iz=0
      iz0=0
      cxymax=0.
      do i=n1,n2,nll 
        iz=iz+1
        if (i.eq.nx1) iz0=iz
        cc=0.
        cy=0.
        cx=0.
        do in=0,npilot-1
          ix=nx1+in
          iy=i+in
c if no overlap goto 55 or 56
          if (iy.lt.1) goto 55
          if (iy.gt.ny) goto 56
c cross-corr:
          yiy=y(iy)
          xix=x(ix)
          cc=cc+x(ix)*yiy
          cy=cy+yiy*yiy
          cx=cx+xix*xix
 55       continue
        enddo
 56     continue
        cxcy=cx*cy
        if (cxcy.gt.cxymax) cxymax=cxcy
        cxxyy(iz)=cx*cy
        ccsig(iz)=cc
      enddo
      nz=iz

      rcxxyy=cxymax*rlim
      nnorm=0
      do iz=1,nz
        cyyi=cxxyy(iz)
        if (cyyi.lt.rcxxyy) then 
          z(iz)=0.0
        else
          rnorm=sqrt(cyyi)
          z(iz)=ccsig(iz)/rnorm
        endif
      enddo

      if (iz0.ne.0) write(*,*)z(iz0),' Zero LAG CCGN'

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE sactraces(name1,name2,sig1,sig2,nlen1,nlen2,
     &  max,dt2,beg2,gcarc2,evdp2,kstnm2,baz2,s1la,s1lo,s2la,s2lo,lnn)
c
c   Returns time series sig1,sig2,dt1,...
c
c--------------------------------------------------------------
c     parameter (max=8192)
c     parameter (max=32768)
c     parameter (max=65536)
c     parameter (max=131072)
c     parameter (max=262144)
c     parameter (max=524288)
c     parameter (max=1048576)
      real sig1(1),sig2(1)
      real tstart1,tend1,gcarc1,evdp1,baz1
      character*80 name1,name2
      integer nsmpl,nlen1,nlen2,max2
      character*8 kstnm2
      logical lnn

      nax=max
      if (nlen2.gt.0) nax=nlen2
      call rsac1(name2,sig2,nlen2,beg2,dt2,nax,nerr)
cc    if (nerr.ne.0) stop 'Error in reading 2nd sac file'
      if (lnn) then
        call getfhv('stla',s2la,nerr)
        call getfhv('stlo',s2lo,nerr)
      endif
      nax=max
      if (nlen1.gt.0) nax=nlen1
      call rsac1(name1,sig1,nlen1,beg1,dt1,nax,nerr)
cc    if (nerr.ne.0) stop 'Error in reading 1st sac file'
c     call getfhv('gcarc',gcarc1,nerr)
c     if (nerr.ne.0) gcarc1=999.
c     call getfhv('evdp',evdp1,nerr)
c     if (nerr.ne.0) evdp1=999.
c     call getkhv('kstnm',kstnm2,nerr)
c     if (nerr.ne.0) kstnm2='XXXX'
c     call getfhv('baz',baz1,nerr)
c     if (nerr.ne.0) baz1=999.
      if (lnn) then
        call getfhv('stla',s1la,nerr)
        call getfhv('stlo',s1lo,nerr)
      endif

c check beg:
      if (abs(beg2-beg1).gt.dt2) then
        write(*,*)' CHECK  sac beg-times: ',beg1, beg2
        stop ' SACTRACES: traces with different beg.'
      endif
c check dt:
      if (dt1.ne.dt2) then
        write(*,*)' CHECK  dt1,dt2: ',dt1,dt2
        stop ' SACTRACES: traces with different dt.'
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE getpar(name1,name2,dt,pow,liasc,lisac,libin,
     &  loasc,losac,lobin,lpcc,lccgn,lenv,lnn,lpow,
     &  nlen1,nlen2,n11,n12,t11,t12,nl1,nl2,tl1,tl2,nll,mmax)
      character*80 name1, name2
      character par*68, par1*68
      logical libin,liasc,lisac
      logical lobin,loasc,losac
      logical lccgn,lpcc,lenv,lnn,lpow
      real pow,dt,tl1,tl2,t11,t12
      integer nlen1,nlen2,n11,n12,n21,n22,nl1,nl2,nll

c  DEFAULTS:
      pow=1.
      dt=1
      lccgn=.false.
      lpcc=.false.
      lenv=.false.
      libin=.false.
      lobin=.false.
      liasc=.false.
      loasc=.false.
      lisac=.true. 
      losac=.false.
      lnn=.false.
      lpow=.false.
      nlen1=0
      nlen2=0
      n11=0
      n12=0
      t11=0.
      t12=0.
      n21=0
      n22=0
      nl1=0
      nl2=0
      tl1=0.
      tl2=0.
      nll=1

cccccccccccccccccccccc
c get the input files:
cccccccccccccccccccccc
c number of arguments: narg
      narg=iargc()
 
      if(narg.lt.3) call usepcc(mmax)

      call getarg(1,name1)
      call getarg(2,name2)

      do iarg=3,narg
        call getarg(iarg,par)
        l=leng(par)
        if (par(1:3).eq.'dt=') then
           par1=par(4:l)
           read(par1,*)dt
        else if (par(1:6).eq.'nlen1=') then
           par1=par(7:l)
           read(par1,*)nlen1
        else if (par(1:6).eq.'nlen2=') then
           par1=par(7:l)
           read(par1,*)nlen2
        else if (par(1:4).eq.'tl1=') then
           par1=par(5:l)
           read(par1,*)tl1
        else if (par(1:4).eq.'tl2=') then
           par1=par(5:l)
           read(par1,*)tl2
        else if (par(1:4).eq.'nl1=') then
           par1=par(5:l)
           read(par1,*)nl1
        else if (par(1:4).eq.'nl2=') then
           par1=par(5:l)
           read(par1,*)nl2
        else if (par(1:4).eq.'nll=') then
           par1=par(5:l)
           read(par1,*)nll
        else if (par(1:4).eq.'t11=') then
           par1=par(5:l)
           read(par1,*)t11
        else if (par(1:4).eq.'t12=') then
           par1=par(5:l)
           read(par1,*)t12
        else if (par(1:4).eq.'n11=') then
           par1=par(5:l)
           read(par1,*)n11
        else if (par(1:4).eq.'n12=') then
           par1=par(5:l)
           read(par1,*)n12
c       else if (par(1:4).eq.'n21=') then
c          par1=par(5:l)
c          read(par1,*)n21
c       else if (par(1:4).eq.'n22=') then
c          par1=par(5:l)
c          read(par1,*)n22
        else if (par(1:4).eq.'pow=') then
           par1=par(5:l)
           read(par1,*)pow
           lpow=.true.
        else if (par(1:4).eq.'pov=') then
           par1=par(5:l)
           read(par1,*)pow
           lpow=.false.
        else if (par(1:4).eq.'ibin') then
           libin=.true.
        else if (par(1:4).eq.'iasc') then
           liasc=.true.
        else if (par(1:4).eq.'isac') then
           lisac=.true.
        else if (par(1:4).eq.'obin') then
           lobin=.true.
        else if (par(1:4).eq.'oasc') then
           loasc=.true.
        else if (par(1:4).eq.'osac') then
           losac=.true.
        else if (par(1:4).eq.'pcc') then
           lpcc=.true.
        else if (par(1:4).eq.'ccgn') then
           lccgn=.true.
        else if (par(1:3).eq.'env') then
           lenv=.true.
        else if (par(1:3).eq.'nn') then
           lnn=.true.
        else
           call usepcc(mmax)
        endif
      enddo

      if (liasc.or.libin) lisac=.false.
      if (.not.lisac.and.nlen1.eq.0) stop 'nlen1=? for none sac data.'
      if (.not.lisac.and.nlen2.eq.0) stop 'nlen2=? for none sac data.'
      if (n11.ne.0.and.t11.ne.0.) stop 'set either n11= or t11='
      if (n12.ne.0.and.t12.ne.0.) stop 'set either n12= or t12='
      if (nl1.ne.0.and.tl1.ne.0.) stop 'set either nl1= or tl1='
      if (nl2.ne.0.and.tl2.ne.0.) stop 'set either nl2= or tl2='

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE envel(x,y,n)
c
c     parameter (nmax=8192)
c     parameter (nmax=32768)
c     parameter (nmax=65536)
c     parameter (nmax=131072)
c     parameter (nmax=262144)
c     parameter (nmax=524288)
c     parameter (nmax=1048576)
      parameter (nmax=8388608)
      real x(nmax),y(nmax)
      complex c(nmax),cc
c
      if (n.gt.nmax) stop ' Change dimensions in envel!'
c
      do i=1,nmax
        c(i)=cmplx(0.,0.)
        y(i)=0.
      enddo
c
      npow=1
      nsmp=2
 5    if (n.le.nsmp) goto 10
      npow=npow+1
      nsmp=nsmp*2
      goto 5
 10   nsmp2=nsmp/2
c
      do 20 i=1,n
 20   c(i)=cmplx(x(i),0.)
        do 30 i=n+1,nsmp
 30   c(i)=cmplx(0.,0.)
c
      dt=1.
      call clogc(npow,c,1.,dt)
c
      do  i=1,nsmp2
        c(i)=2.*c(i)
        c(nsmp-i+2)=cmplx(0.,0.)
      enddo
      c(1)=cmplx(0.,0.)
      c(nsmp2+1)=cmplx(0.,0.)
c
      call clogc(npow,c,-1.,dt)
c
      do i=1,n
        cc=c(i)
        y(i)=cabs(cc)
      enddo
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE wrsac(name1,kstnm,b,dt,gcarc,evdp,y,nsmpl,
     & lnn,s1la,s1lo,s2la,s2lo)

      real y(1)
      character*8 kstnm
      character*(*) name1
      logical lnn

      call newhdr
      call setnhv('npts',nsmpl,nerr)
      call setfhv('delta',dt,nerr)
      call setkhv('kstnm',kstnm,nerr)
      if (lnn) then
        call setfhv('evla',s1la,nerr)
        call setfhv('evlo',s1lo,nerr)
        call setfhv('stla',s2la,nerr)
        call setfhv('stlo',s2lo,nerr)
      else
        call setfhv('gcarc',gcarc,nerr)
        call setfhv('evdp',evdp,nerr)
c       call setfhv('evla',evla,nerr)
c       call setfhv('evlo',evlo,nerr)
c       call setfhv('stla',stla,nerr)
c       call setfhv('stlo',stlo,nerr)
      endif

      e=b+(nsmpl-1)*dt
      o=0.

c     call setfhv('e',e,nerr)
c     call setfhv('o',o,nerr)
      call setfhv('b',b,nerr)
      call setihv('iztype','io',nerr)
c
      call wsac0(name1,dum,y,nerr)
      if(nerr.ne.0) stop 'Error in writing output file'
      
      return
      end
c__________________________________________________________________
      SUBROUTINE usepcc(mmax)
      integer mmax
      write(*,*)' ' 
      write(*,*)'Phase Cross-Correlations:'
      write(*,*)'========================='
      write(*,*)' '
      write(*,*)'USAGE: pcc5g in1 in2 parameters'
      write(*,*)'  '
      write(*,*)'PARAMETERS:'
      write(*,*)' nl1= nl2= n11= n12= nlen1= nlen2='  
      write(*,*)' pow= iasc ibin isac oasc obin osac ccgn env nn'
      write(*,*)' '
      write(*,*)' in1 in2       : input data (pilot & trace)' 
      write(*,*)' iasc,ibin,isac: input format of in1,in2'
      write(*,*)' nlen1,nlen2   : # of samples' 
      write(*,*)'                 (if not set: use entire t-series.)'
      write(*,*)' nl1,nl2   : relative sample lag.'
      write(*,*)' tl1,tl2   : relative time lag in seconds.'
      write(*,*)'             (use either nl1=,nl2= or tl1=,tl2=)'
      write(*,*)' nll       : sample lag interval (nll*dt, nll>=1)'
      write(*,*)'             Only useful if dt much smaller than '
      write(*,*)'             accuracy or for quick and dirty tests.'
      write(*,*)' n11,n12   : pilot start/end samples (1st t-series).'
      write(*,*)' t11,t12   : pilot start/end time in secondes'//
     & ' (1st t-series)'
      write(*,*)'             (use either n11=,n12= or t11=,t12=).'
      write(*,*)' pow       : pcc power (either pow or pov)'
      write(*,*)'             ... (Sum)*|Sum|**(pow-1)'
      write(*,*)' pov       : pcc power (either pow or pov)'
      write(*,*)'             ... Sum_of( |a+b|**pov -|a-b|**pov )'
      write(*,*)' ccgn      : geom. norm. cross-corr. is computed.'
      write(*,*)' pcc       : phase cross-corr is computed.'
      write(*,*)' oasc,obin,osac: output format (any combination)'
      write(*,*)' env       : output also envelopes'
      write(*,*)' nn        : set evla/evlo stla/stlo on sac output '
      write(*,*)'             using stla/stlo from in1 & in2, '//
     & 'respectively.'
      write(*,*)'             evla=stla1 evlo=stlo1 stla=stla2 '//
     & 'stlo=stlo2'
      write(*,*)'             I.e., SAC computes interstation '//
     & 'distance.'
      write(*,*)' '
      write(*,*)'DEFAULTS: '
      write(*,*)' nl1=0 nl2=0 nll=1 isac pov=1 n11=1 n12=nlen1 '
      write(*,*)' '
      write(*,*)'DIMENSIONS: '
      write(*,*)' maximum length of t series (samples):',mmax
      write(*,*)' '
      write(*,*)'EXAMPLE1: '
      write(*,*)'(pcc, sac output, t-lag:-20s to 20s)'
      write(*,*)' pcc5g STA1.sac STA2.sac pcc tl1=-20 '//
     & 'tl2=20 osac'
      write(*,*)'EXAMPLE2: '
      write(*,*)' (pcc,ccgn, lag:-200 to 200 samples, sac output,'//
     & ' write stla/stlo'
      write(*,*)' to SAC-header, pcc power pov=2)'
      write(*,*)' pcc5g STA1.sac STA2.sac pcc ccgn nl1=-200 '//
     & 'tl2=200 osac nn pov=2'
      write(*,*)'EXAMPLE3: '
      write(*,*)' (like EXAMPLE 2, but use only STA1.sac, '//
     & 'pilot from 10s to 123.2s)'
      write(*,*)' pcc5g STA1.sac STA1.sac pcc ccgn nl1=-200 '//
     & 'tl2=200 osac nn t11=10 t12=123.2'
      write(*,*)
      write(*,*)'IMPORTANT: '
      write(*,*)' i)   t-series must have same start time.'
      write(*,*)' ii)  A pilot can be selected out of the 1st trace.'
      write(*,*)' iii) Time series can have same length.'
      write(*,*)' iv)  Correlation lag interval equals the sampling'
      write(*,*)'      interval. If accuracy needs to be increased then'
      write(*,*)'      interpolate the t-series to a smaller sample'
      write(*,*)'      interval dt (e.g., in SAC).'
      write(*,*)' v)   Amplitudes of t-series should not bt too small.'
      write(*,*)'      If amplitudes become < 10-07 then just multiply'
      write(*,*)'      your data with a constant factor.'
      write(*,*)' '
      write(*,*)'LAST CHANGE: 09/12/2014'
      write(*,*)'Questions,bugs,improvements,... => '//
     & 'schimmel@ictja.csic.es'
      write(*,*)' '
      stop
      return
      end
c__________________________________________________________________
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
c--- at t=0, the spectral density is dt for all values of the frequency.
c
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

c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------
      SUBROUTINE analytic(s1,c1,nsamp)
c
c Compute analytic signal c1(nsamp) from real time series s1(nsamp).
c Frequency domain computation!
c
      real s1(1)
      complex c1(1)
c
      dt=1.
c
c Determine power for FFT:
      npow=1
      nsmp=2
 5    if (nsamp.le.nsmp) goto 10
      npow=npow+1
      nsmp=nsmp*2
      goto 5
 10   nsmp2=nsmp/2
c
      do 20 i=1,nsamp
 20   c1(i)=cmplx(s1(i),0.)
        do 30 i=nsamp+1,nsmp
 30   c1(i)=cmplx(0.,0.)
c
c FFT:
      call clogc(npow,c1,1.,dt)
c
c calculating the analytical signal:
      do i=2,nsmp2
          c1(i)=2.*c1(i)
          c1(nsmp+2-i)=cmplx(0.,0.)
      enddo
c       c1(1)=cmplx(0.,0.)
c     c1(nsmp+1)=cmplx(0.,0.)
c
c IFFT:
      call clogc(npow,c1,-1.,dt)
c
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE pcc_time_exp_out3(s1,r1,r2,nlag1,nlag2,nll,wu,
     &  n1,n2,ll,n11)
c
c     TIME DOMAIN PCC (BSSA,1999). Function returns real value as 
c     function of relative sample lag (see below).
c
c     r1 & r2 contain the phasors (envelope normalized analytic signals) 
c     for the pilot and trace.
c     n1 & n2 are the corresponding number of samples.
c     r1 is shifted with respect to r2, i.e., a positive lag means
c     that r1 is shifted to larger time. The time of r1(1) equals
c     the time at r2(n11) and nlag=n11+nl, where nl is the sample
c     lag (positive or negative integer value). nl=0 means zero lag.
c      
c     Partial overlap is permitted and depends on nlag.
c     Complete overlap for 0 <= nlag < (n2-n1+1) : 
c     pcc[ r1(i) , r2(i+nlag) ] for i=1,..,i21
c     wu is the power of PCC. (wu=1 corresponds to PCC from
c     Schimmel, BSSA, 1999.)
c
c     last change: 14/12/2011 schimmel@ictja.csic.es
c     last change: 28/05/2013 schimmel@ictja.csic.es
c
c     real wu,r1(1),r2(1),cosd,rdum,difarg
c     integer n1,n2,nnorm,nlag,i,nn

      complex r1(1),r2(1),c1,c2
      real s1(1)

      ll=0
      do nlag=nlag1,nlag2,nll

      ll=ll+1
      if (nlag.lt.1-n1)  stop 'No overlap at all!'
      if (nlag.ge.n2) stop 'No overlap at all!'

c 1st sample on 1st trace:
      i11=2-nlag
      if (i11.lt.1) i11=1
c 1st sample on 2nd trace:
      i12=nlag
      if (i12.lt.1) i12=1
c last sample on 1st trace 
      nn2=n2-i12
      nn1=n1-i11
      i21=min(nn1,nn2)+i11
      
      nn=i12
      d1=0.
      d2=0.
      none=0
      do i=i11,i21
        c1=r2(nn)
        c2=r1(i)
        if (cabs(c1)*cabs(c2).ne.0) then
          d1=d1+cabs(c1+c2)
          d2=d2+cabs(c1-c2)
        else
          none=none+1
        endif
        nn=nn+1
      enddo
      rdum=d1-d2

      nnorm=i21-i11+1-none
      rdorm=rdum*.5/float(nnorm)
      if (wu.eq.1) then
         rdum=rdorm
      elseif (wu.eq.2) then
         rdum=rdorm*abs(rdorm)
      else
         rdum=rdorm*abs(rdorm)**(wu-1)
      endif

      s1(ll)=rdum

c set outside loop.
      if (nlag.eq.n11) write(*,*)rdum,' Zero LAG PCC'

      enddo

      return
      end
c____________________________________________________________________
      SUBROUTINE pcc_time_exp_in3(s1,r1,r2,nlag1,nlag2,nll,wu,
     &  n1,n2,ll,n11)
c
c     TIME DOMAIN PCC (BSSA,1999). Function returns real value as 
c     function of relative sample lag (see below).
c
c     r1 & r2 contain the phasors (envelope normalized analytic signals) 
c     for the pilot and trace.
c     n1 & n2 are the corresponding number of samples.
c     r1 is shifted with respect to r2, i.e., a positive lag means
c     that r1 is shifted to larger time. The time of r1(1) equals
c     the time at r2(n11) and nlag=n11+nl, where nl is the sample
c     lag (positive or negative integer value). nl=0 means zero lag.
c      
c     Partial overlap is permitted and depends on nlag.
c     Complete overlap for 0 <= nlag < (n2-n1+1) : 
c     pcc[ r1(i) , r2(i+nlag) ] for i=1,..,i21
c     wu is the power of PCC. (wu=1 corresponds to PCC from
c     Schimmel, BSSA, 1999.)
c
c     last change: 14/12/2011 schimmel@ictja.csic.es
c     last change: 28/05/2013 schimmel@ictja.csic.es
c
c     real wu,r1(1),r2(1),cosd,rdum,difarg
c     integer n1,n2,nnorm,nlag,i,nn

      complex r1(1),r2(1),c1,c2,cc1,cc2
      real s1(1)

      rr=.5**wu
      ll=0
      do nlag=nlag1,nlag2,nll

      ll=ll+1
      if (nlag.lt.1-n1)  stop 'No overlap at all!'
      if (nlag.ge.n2) stop 'No overlap at all!'

c 1st sample on 1st trace:
      i11=2-nlag
      if (i11.lt.1) i11=1
c 1st sample on 2nd trace:
      i12=nlag
      if (i12.lt.1) i12=1
c last sample on 1st trace 
      nn2=n2-i12
      nn1=n1-i11
      i21=min(nn1,nn2)+i11

      nn=i12
      rdum=0.
      d1=0.
      d2=0.
      
      none=0
      if (wu.eq.1) then
        do i=i11,i21
          c1=r2(nn)
          c2=r1(i)
          if (cabs(c1)*cabs(c2).ne.0.) then
            d1=d1+cabs(c1+c2)
            d2=d2+cabs(c1-c2)
          else
            none=none+1
          endif
          nn=nn+1
        enddo
        rdum=d1-d2
        nnorm=i21-i11+1-none
        rnorm=0.5/float(nnorm)
        rdum=rdum*rnorm
      else if (wu.eq.2) then
        do i=i11,i21
          rrc1=real(r2(nn))
          ric1=aimag(r2(nn))
          rrc2=real(r1(i))
          ric2=aimag(r1(i))
          if ((rrc1+ric1)*(rrc1+ric1).ne.0) then
            rdum=rdum+rrc1*rrc2+ric1*ric2
          else
            none=none+1
          endif
          nn=nn+1
        enddo
        nnorm=i21-i11+1-none
        rdum=rdum/float(nnorm)
      else
        do i=i11,i21
          c1=r2(nn)
          c2=r1(i)
          if (cabs(c1)*cabs(c2).ne.0.) then
            d1=cabs(c1+c2)
            d2=cabs(c1-c2)
            rdum=rdum+(d1**wu-d2**wu)
          else
            none=none+1
          endif
          nn=nn+1
        enddo
        nnorm=i21-i11+1-none
        rnorm=rr/float(nnorm)
        rdum=rdum*rnorm
      endif
    

      s1(ll)=rdum

c set outside loop.
      if (nlag.eq.n11) write(*,*)rdum,' Zero LAG PCC'

      enddo

      return
      end
c____________________________________________________________________
