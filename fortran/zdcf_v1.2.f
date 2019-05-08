ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c       The Z-transformed Discrete Correlation Function algorithm (ZDCF)      c
c                                                                             c
c                                Version 1.2                                  c
c                                                                             c
c   Reference:                                                                c
c                                                                             c
c   T. Alexander, 1997, 'Is AGN variability correlated with other AGN         c
c   properties? - ZDCF analysis of small samples of sparse light curves',     c
c   in "Astronomical Time Series", eds. D. Maoz, A. Sternberg and             c
c   E.M. Leibowitz, Kluwer, Dordrecht, p. 163                                 c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c SAMPLE RUN (user input marked by "<-USER". See explanations below)          c
c ----------                                                                  c
c                                                                             c
c ZDCF V1.2 begins.                                                           c
c Auto-correlation or cross-correlation? (1/2):                               c
c 1     <-USER                                                                c
c Enter output files prefix:                                                  c
c auto_g <-USER                                                               c
c Uniform sampling of light curve? (y/n):                                     c
c n     <-USER                                                                c
c Enter minimal number of points per bin (0 for default):                     c
c 0     <-USER                                                                c
c Omit zero-lag points? (y/n):                                                c
c y     <-USER                                                                c
c How many Monte Carlo runs for error estimation?                             c
c 100   <-USER                                                                c
c Enter name of 1st light curve file:                                         c
c g     <-USER                                                                c
c Enter name of 2nd light curve file:                                         c
c g_2   <-USER                                                                c
c                                                                             c
c dbg.lc1  written (contains   59 points) ...                                 c
c dbg.lc2  written (contains   63 points) ...                                 c
c                                                                             c
c =========================================================================== c
c ZDCF PARAMETERS:                                                            c
c Autocorrelation?    F                                                       c
c Uniform sampling?   F                                                       c
c Omit zero lags?     T                                                       c
c Minimal # in bin:   11                                                      c
c # of MonteCarlo:    100                                                     c
c MonteCarlo seed:    123                                                     c
c   59 points in inp1                                                         c
c   63 points in inp2                                                         c
c =========================================================================== c
c                                                                             c
c Binning with minimum of  11 points per bin and resolution of 1.00E-03 .     c
c                                                                             c
c   71 bins actually used,   307 inter-dependent pairs discarded.             c
c                                                                             c
c   tau       -sig(tau)  +sig(tau)   dcf        -err(dcf)  +err(dcf) (#bin)   c
c                                                                             c
c  -5.346E+01  2.233E+01  4.374E-01  4.480E-01  3.394E-01  2.459E-01 (  10)   c
c  -5.104E+01  9.915E-01  1.058E-01  6.378E-01  2.524E-01  1.636E-01 (  11)   c
c   .                                                      .                  c
c   .                                                      .                  c
c   .                                                      .                  c
c   5.203E+01  1.396E-01  9.683E-01  5.681E-02  3.433E-01  3.304E-01 (  11)   c
c   5.498E+01  1.020E+00  2.047E+01  1.553E-01  3.704E-01  3.313E-01 (  10)   c
c                                                                             c
c dbg.dcf written...                                                          c
c Program ended.                                                              c
c                                                                             c
c EXPLANATION                                                                 c
c -----------                                                                 c
c                                                                             c
c o The input light curve files should be in 3 columns (free, floating        c
c   point format):                                                            c
c   time (ordered); flux / magnitude; absolute error on flux / magnitude.     c
c                                                                             c
c o The sign of the time lag is defined as                                    c
c   tau = t(2nd light curve) - t(1st light curve)                             c
c                                                                             c
c o The program creates 3 output files with the same prefix. dbg.lc1, dbg.lc2 c
c   are the 2 light curves (in the same format as the input) after points     c
c   with identical times are averaged. dbg.dcf is the resulting ZDCF in 7     c
c   columns, as in the example above (note: the extreme lags may be based     c
c   on less than the minimal 11 points. these are the "left-over" points and  c
c   are used just for giving an "impression" of the trends at extreme lags.   c
c   However, these ZDCF points are NOT statistically reliable).               c
c                                                                             c
c o If the light curve sampling IS uniform, it is possible to force the       c
c   program not to collect different time lags into a single bin              c
c                                                                             c
c o The default minimal number of points per bin is 11.                       c
c                                                                             c
c o It is recommended to omit the zero-lag points.                            c
c                                                                             c
c o 100 Monte Carlo runs seem to be sufficient in most cases.                 c
c                                                                             c
c FORTRAN ISSUES                                                              c
c --------------                                                              c
c                                                                             c
c o This implementation is designed for speed at the cost of large memory     c
c   requirements. The parameters MPNTS and MBINS in the main's variable       c
c   declaration header can be changed to fit your system's limitations.       c
c   (Here and below, Fortran words are emphasized in upper-case, although     c
c   they appear in the program in lower-case).                                c
c                                                                             c
c o Non-standard dynamic allocation functions, MALLOC and FREE, are used in   c
c   the subroutine ALCBIN. If your system does not support MALLOC, un-comment c
c   the declarations currently bracketed by C-NO-MALLOC and comment the       c
c   declarations, MALLOC calls and FREE calls which are currently bracketed   c
c   by C-MALLOC comments.                                                     c
c                                                                             c
c o The small parameter epsilon, which controls the binning of pairs with     c
c   close lags (see section 2.2.2 in the paper), is currently not user-       c
c   adjustable. It can be changed by modifying the parameter RSLUTN in the    c
c   subroutine ALCBIN. Epsilon = RSLUTN * (The maximal time lag in the data). c
c   Currently RSLUTN is set to 0.001 .                                        c
c                                                                             c
c o The random number generation routine rndnrm assumes the existence of a    c
c   function call x = rand(0), which produces a pseudo random real number x   c
c   in the interval [0..1], and is initialized by x = rand(iseed), where      c
c   iseed is an integer > 1.                                                  c
c                                                                             c
c o The seed for the Monte Carlo random simulations, SEED, is initialized by  c
c   a DATA statement at the end of the main's variable declaration header.    c
c                                                                             c
c VERSION LOG                                                                 c
c -----------                                                                 c
c                                                                             c
c V1.0 ??/??/95: Original version                                             c
c V1.1 09/05/96: Minor corrections for compatibility with Alpha/OSF           c
c V1.2 21/10/99: Use open source routines for random numbers and sorting      c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer mpnts,mbins,minbin
c*****************************************************************************c
c                                                                             c
c     SETTING THE MAXIMAL NO. OF POINTS IN THE LIGHT CURVES (mpnts)           c
c     AND THE MAXIMAL NO. OF BINS IN THE ACF / CCF (mbins).                   c
c     ATTENTION: when not using alloc, mpnts should be also defined with      c
c                the same value in the function alcbin.                       c
c                                                                             c
      parameter (mpnts = 1550, mbins = 500)
c*****************************************************************************c
      parameter (minbin=mpnts)
      integer i,n
      integer na,nb,nbins,used,unused
      integer minpts
      logical autocf,no0lag,nwdata,unsmpl
      real ta(mpnts),tb(mpnts)
      real a(mpnts),b(mpnts)
      real erra(mpnts),errb(mpnts)
      real tdcf(mbins),sigtm(mbins),sigtp(mbins),dcf(mbins)
      real sigdcm(mbins),sigdcp(mbins)
      logical wuseda(mpnts),wusedb(mpnts)
      integer inbin(mbins),enough
      parameter (enough=11)
c Areas for Monte Carlo estimates of the measurement errors
      integer ncarlo
      real acarlo(mpnts),bcarlo(mpnts)
      real avz(mbins)
c Work areas for clcdcf
      real wa(minbin),wb(minbin)
      real wta(minbin),wtb(minbin)
      real werra(minbin),werrb(minbin)
      integer wai(minbin,mbins),wbj(minbin,mbins)
      real vtau(minbin,mbins)
c     
      character*72 infil1,infil2,outfil,prefix
      integer oerr,cerr,lprfix
      integer seed
      real expz,sigz,eps
      parameter (eps=1e-7)
c     
      integer inpi
      character*1 inpc
      real fishe,fishs
      real rand
      real z,r
      data seed/123/     
c     
c Reading program parameters
c     
      print *,'ZDCF V1.2 begins.'
      print *,'Auto-correlation or cross-correlation? (1/2):'
c     read *,inpi
      inpi = 1
      print *,inpi
      autocf = (inpi .eq. 1)
c     
      print *,'Enter output files prefix:'
      read '(a72)',prefix
      do i = 1,len(prefix)
         if (prefix(i:i) .eq. ' ') then
            lprfix = i-1
            goto 200
         endif
      enddo
200   continue
c     
      print *,'Uniform sampling of light curve? (y/n):'
c     read '(a1)',inpc
      inpc = 'n'
      print *,inpc
      unsmpl = (inpc .eq. 'y' .or. inpc .eq. 'Y')
      print *,'Enter minimal number of points per bin (0 for default):'
c     read *,minpts
      minpts = 11
      print *,minpts
c if minpts = 0, use default
      if (minpts .le. 0) minpts = enough
c     
      print *,'Omit zero-lag points? (y/n):'
c     read '(a1)',inpc
      inpc = 'y'
      print *,inpc
      no0lag = (inpc .eq. 'y' .or. inpc .eq. 'Y') 
c     
      print *,'How many Monte Carlo runs for error estimation?'
c     read *,ncarlo
      ncarlo = 100
      print *,ncarlo
      if (ncarlo.le.1) ncarlo = 0
c     
c Reading the observed data
c     
      print *,'Enter name of 1st light curve file:'
c     read '(a72)',infil1
      infil1 = trim(prefix) // '.csv'
      print *,infil1
      if (.not. autocf) then
         print *,'Enter name of 2nd light curve file:'
c        read '(a72)',infil2
         infil2 = trim(prefix) // '_2.csv'
      endif
      call redobs(na,nb,ta,tb,a,b,erra,errb,mpnts,
     >            autocf,infil1,infil2)
c     
c Writing the condensed light curve data 
c     
      print *
      outfil = prefix(:lprfix)//'.lc1'
      open (unit=11,file=outfil,status='UNKNOWN',iostat=oerr,err=901)
      write (11,'(1p,3(1x,e12.5))')
     >   (ta(i),a(i),erra(i),i=1,na)
      close (unit=11,status='KEEP',iostat=cerr,err=902)
      print *,outfil(:lprfix+5),' written (contains ',na,' points) ...'
c     
      outfil = prefix(:lprfix)//'.lc2'
      open (unit=12,file=outfil,status='UNKNOWN',iostat=oerr,err=901)
      write (12,'(1p,3(1x,e12.5))')
     >   (tb(i),b(i),errb(i),i=1,nb)
      close (unit=12,status='KEEP',iostat=cerr,err=902)
      print *,outfil(:lprfix+5),' written (contains ',nb,' points) ...'
c     
      print *
      print *,('=',i=1,75)
      print *,'ZDCF PARAMETERS:'
      print *,'Autocorrelation?  ',autocf
      print *,'Uniform sampling? ',unsmpl
      print *,'Omit zero lags?   ',no0lag
      print *,'Minimal # in bin: ',minpts
      print *,'# of MonteCarlo:  ',ncarlo 
      print *,'MonteCarlo seed:  ',seed 
      print *,na,' points in ',infil1(1:40)
      if (.not. autocf) print *,nb,' points in ',infil2(1:40)
      print *,('=',i=1,75)
      print *
c     
c Estimating the effects of the measurement errors by Monte Carlo simulations
c     
      i = rand(seed)
      if (ncarlo .gt. 1) then
         do i = 1,mbins
            avz(i) = 0.0
         enddo
         nwdata = .true.
         do i = 1,ncarlo
            call simerr(a,acarlo,erra,na)
            if (autocf) then
               call clcdcf(ta,acarlo,erra,na,
c                          i  i      i    i
     >                     tb,acarlo,errb,nb,
c                          i  i      i    i
     >             autocf,no0lag,unsmpl,used,unused,nwdata,
c                  i      i      i      o    o      i/o    
     >             minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp,
c                  i/o    o    o     o     o   o      o    
     >             inbin,minbin,nbins,mbins,
c                  i    i      i/o   i
     >             wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb,
c                  w   w   w  w  w     w     w   w   w      w
     >             vtau)
c                  o
            else
               call simerr(b,bcarlo,errb,nb)
c     
c calculating the discrete correlation function for Monte Carlo error approx.
c     
               call clcdcf(ta,acarlo,erra,na,
c                          i  i i    i
     >                     tb,bcarlo,errb,nb,
c                          i  i i    i
     >             autocf,no0lag,unsmpl,used,unused,nwdata,
c                  i      i      i      o    o      i/o    
     >             minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp,
c                  i/o    o    o     o     o   o      o      
     >             inbin,minbin,nbins,mbins,
c                  i    i      i/o   i     
     >             wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb,
c                  w   w   w  w  w     w     w   w   w      w
     >             vtau)
c                  o
            endif
            if (i .eq. 1) then
               print *
               print *,nbins,' bins actually used, ',unused,
     >              ' inter-dependent pairs discarded.'
            endif
c     
c The summing and averaging is done in z-space.
c     
            do n = 1,nbins
               r = dcf(n)
               if (r .lt. -1.0+eps) then
                  r = -1.0+eps
               else if (r .gt.  1.0-eps) then
                  r =  1.0-eps
               endif   
               avz(n) = avz(n)+log((1.+r)/(1.-r))/2.
            enddo
         enddo
         do i = 1,nbins
            z = avz(i)/(ncarlo+0.)
            n = inbin(i)
            dcf(i) = tanh(z) 
            r = dcf(i)
            if (r .lt. -1.0+eps) then
               r = -1.0+eps
            else if (r .gt. 1.0-eps) then
               r =  1.0-eps
            endif
            z = log((1.+r)/(1.-r))/2.0
            sigz = fishs(r,n)
            expz = fishe(r,n)
            sigdcm(i) = r-tanh(expz-sigz)
            sigdcp(i) = tanh(expz+sigz)-r
         enddo
      else
c     
c calculating the discrete correlation function w/o Monte Carlo error approx.
c     
         nwdata = .true.
         call clcdcf(ta,a,erra,na,
c                    i  i i    i
     >               tb,b,errb,nb,
c                    i  i i    i
     >               autocf,no0lag,unsmpl,used,unused,nwdata,
c                    i      i      i      o    o      i/o   
     >               minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp,
c                    i/o    o    o     o     o   o      o      
     >               inbin,minbin,nbins,mbins,
c                    o     i      i/o   i,   
     >          wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb,
c               w   w   w  w  w     w     w   w   w      w
     >          vtau)
c               o
         print *
         print *,nbins,' bins actually used, ',unused,
     >           ' inter-dependent pairs discarded.'
      endif
c     
c printing the results: tdcf, dcf, sigdcf
c     
      outfil = prefix(:lprfix)//'.dcf'
      open (unit=13,file=outfil,status='UNKNOWN',iostat=oerr,err=901)
      print *
      print '(1p,6(1x,a10),'' ('',a4,'')'')',
     >      ' tau      ','-sig(tau) ','+sig(tau) ',' dcf      ',
     >      ' -err(dcf)',' +err(dcf)','#bin'
      print *
      do i = 1,nbins
         print '(1p,6(1x,e10.3),'' ('',i4,'')'')',
     >         tdcf(i),sigtm(i),sigtp(i),
     >         dcf(i),sigdcm(i),sigdcp(i),inbin(i)
         write (13,'(1p,6(1x,e10.3),1x,i4)') 
     >         tdcf(i),sigtm(i),sigtp(i),
     >         dcf(i),sigdcm(i),sigdcp(i),inbin(i)
      enddo
      close (unit=13,status='KEEP',iostat=cerr,err=902)
      print *
      print *,outfil(:lprfix+4),' written...'
c     
      print *,'ZDCF ended.'
      stop
 901  print *,'ERROR ON OPENING FILE - ERRCOD=',oerr
      stop
 902  print *,'ERROR ON CLOSING FILE - ERRCOD=',cerr
      stop
      end      
c//////////////////////////////////////////////////////////////////////////////
c     
c Reading the observed points
c     
      subroutine redobs(na,nb,ta,tb,a,b,erra,errb,mpnts,
     >                  autocf,infil1,infil2)
      integer na,nb
      real ta(mpnts),tb(mpnts),a(mpnts),b(mpnts)
      real erra(mpnts),errb(mpnts)
      logical autocf
      character*72 infil1,infil2
      integer nsamet,oerr,cerr
c     
c Reading the 1st light curve
c     
      open (unit=1,file=infil1,status='UNKNOWN',iostat=oerr,err=901)
      nsamet = 1
      na = 1
 1    continue
         if (na.gt.mpnts) then
            print *,
     >'Only 1st ',mpnts,' points will be read from ',infil1(1:40)
            goto 2
         endif
         read (1,*,end=2) ta(na),a(na),erra(na)
c     
c Averaging observations with identical times
c     
         if (na.gt.1 .and. ta(na).eq.ta(na-1)) then
            a(na-1) = a(na-1)+a(na)
            erra(na-1) = erra(na-1)+erra(na)
            nsamet = nsamet+1
         else if (na.gt.1) then
            a(na-1) = a(na-1)/nsamet
            erra(na-1) = erra(na-1)/nsamet
            nsamet = 1
            na = na+1
         else
            na = na+1
         endif
         goto 1
 2    continue
      close (unit=1,status='KEEP',iostat=cerr,err=902)
      na = na-1
      if (na .le. 0) then
         print *,'No data in ',infil1
         stop
      endif
      a(na) = a(na)/nsamet
      erra(na) = erra(na)/nsamet
c     
c Reading the 2nd light curve
c     
      if (autocf) then
         nb = na
         do i = 1,nb
            tb(i) = ta(i)
            b(i) = a(i)
            errb(i) = erra(i)
         enddo
      else
         open (unit=2,file=infil2,status='UNKNOWN',
     >         iostat=oerr,err=901)
         nsamet = 1
         nb = 1
 3       continue
            if (nb.gt.mpnts) then
               print *,
     >'Only 1st ',mpnts,' points will be read from ',infil2(:40)
               goto 4
            endif
            read (2,*,end=4) tb(nb),b(nb),errb(nb)
c     
c Averaging observations with identical times
c     
            if (nb.gt.1 .and. tb(nb).eq.tb(nb-1)) then
               b(nb-1) = b(nb-1)+b(nb)
               errb(nb-1) = errb(nb-1)+errb(nb)
               nsamet = nsamet+1
            else if (nb.gt.1) then
               b(nb-1) = b(nb-1)/nsamet
               errb(nb-1) = errb(nb-1)/nsamet
               nsamet = 1
               nb = nb+1
            else
               nb = nb+1
            endif
            goto 3
 4       continue
         close (unit=2,status='KEEP',iostat=cerr,err=902)
         nb = nb-1
         if (nb .le. 0) then
            print *,'No data in ',infil2
            stop
         endif
         b(nb) = b(nb)/nsamet  
         errb(nb) = errb(nb)/nsamet
      endif
      return
 901  print *,'ERROR ON OPENING FILE - ERRCOD=',oerr
      return
 902  print *,'ERROR ON CLOSING FILE - ERRCOD=',cerr
      return
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c Generating a "true" signal by substracting gaussian noise from observations.
c The error erra is the ABSOLUTE error in a.
c     
      subroutine simerr(a,acarlo,erra,na)
      implicit none
      integer na
      real a(na),acarlo(na),erra(na)
      integer i
      real rndnrm
c     
      do i = 1,na
         acarlo(i) = a(i)-erra(i)*rndnrm()
      enddo      
      return
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c Allocating points to bins
c     
      subroutine alcbin(ta,na,tb,nb,wuseda,wusedb,
     >                  no0lag,autocf,minpts,
     >                  tdcf,sigtm,sigtp,
     >                  wai,wbj,inbin,nbins,mbins,minbin,unsmpl,
     >                  vtau)
      implicit none
c*****************************************************************************c
      real rslutn
      parameter (rslutn=0.001)
c*****************************************************************************c
      integer na,nb,nbins,mbins,minbin,minpts
      real ta(na),tb(nb)
      integer inbin(mbins)
      real tdcf(mbins),sigtm(mbins),sigtp(mbins)
      integer wai(minbin,mbins),wbj(minbin,mbins)
      real vtau(minbin,mbins)
      logical no0lag,autocf
      logical wuseda(na),wusedb(nb)
c-malloc(
c      integer malloc
c      pointer (iwaidx,waidx)
c      pointer (iwbidx,wbidx)
c      pointer (iwtau,wtau)
c      pointer (iidx,idx)
c      integer idx(1),waidx(1),wbidx(1)
c      real wtau(1)
c-malloc)
c-no-malloc(
      integer mpnts
c     ATTENTION: mpnts here should equal the value defined in the MAIN!!!
      parameter(mpnts=1550)
      integer idx(mpnts*mpnts)
      integer waidx(mpnts*mpnts),wbidx(mpnts*mpnts)
      real wtau(mpnts*mpnts) 
c-no-malloc)
      integer i,j,p,np,pfr,pmax,incr,inb,nnegtv,plo,phi
      real tij,tolrnc
      logical frstym,unsmpl
      integer min0,max0
      data frstym/.true./
c     
c Allocating the dynamic work areas
c     
c-malloc(
c      iwaidx = malloc(na*nb*4)
c      iwbidx = malloc(na*nb*4)
c      iwtau  = malloc(na*nb*4)
c      iidx   = malloc(na*nb*4)
c-malloc)
c     
c Calculate all the time lag points
c     
      if (frstym) then
         frstym = .false.
         print 
     >      '(''Binning with minimum of '',i3,'//
     >      ''' points per bin and resolution of '',1p,g8.2,'' .'')',
     >      minpts,rslutn
      endif
      np = 0
      if (autocf) then
         do i = 1,na
            do j = i,nb
               tij = tb(j)-ta(i)
               if (no0lag .and. tij .eq. 0.0) goto 11
               np = np+1
               wtau(np) = tij
               waidx(np) = i
               wbidx(np) = j
11             continue
            enddo
         enddo
      else
         do i = 1,na
            do j = 1,nb
               tij = tb(j)-ta(i)
               if (no0lag .and. tij .eq. 0.0) goto 12
               np = np+1
               wtau(np) = tij
               waidx(np) = i
               wbidx(np) = j
12             continue
            enddo
 2       enddo
      endif
c     
c Sort according to increasing time lag
c     
      call srtind(wtau,np,idx,1)
c     
c calculating the tolerance level for lags to be considered the same
c     
      tij = wtau(idx(np))-wtau(idx(1))
      tolrnc = tij*rslutn
c     
c Looping on bins
c     
      nbins = 0
c     
c If binned CCF: binning from median time-lag upwards and backwards!
c     
      if (autocf .or. unsmpl) then
         pfr = 1
         pmax = np
         incr = 1
      else
         pfr = np/2
         pmax = 1
         incr = -1
      endif
20    continue
      inb = 0
      tij = incr*1e30
      nbins = nbins+1
      if (nbins .gt. mbins) then
         print *,'Not enough work area - increase mbins in MAIN!'
         stop
      endif
      tdcf(nbins) = 0.0
c     
c Initialize used point flag vectors
c     
      do i = 1,na
         wuseda(i) = .false.
      enddo
      do i = 1,nb
         wusedb(i) = .false.
      enddo
c     
c Collect points into bins that contain at least "enough" points,
c but do not break points with the same lag (up to the tolerance)
c into separate bins
c     
      do i = pfr,pmax,incr
         p = idx(i)
c Check whether bin is full
         if ( ( incr*(wtau(p)-tij-incr*tolrnc) .gt. 0.0 .and. 
     >          inb .ge. minpts ) .or.
     >       (unsmpl .and. incr*(wtau(p)-tij-incr*tolrnc) .gt. 0.0) .or. 
     >       (i .eq. pmax) ) then
c     
c Bin is full: Calculating tau and its std (before proceeding to the next 
c bin, at label 20)  
c     
            inbin(nbins) = inb
            tdcf(nbins) = tdcf(nbins)/inb
            if (unsmpl) then
               sigtm(nbins) = 0.0
               sigtp(nbins) = 0.0
               pfr = i
c If not enough points in bin, ignore it
               if (inb .lt. minpts) nbins = nbins-1
               if (pfr .ne. pmax) goto 20
            else
c If the last point is alone in its bin, we ignore it (to avoid subsequent
c divisions by zero) and get immediately out of loop (label 40)
               if (inb .le. 1) then
                  nbins = nbins-1
                  goto 40
               endif
c     
c Finding the 0.3413 (+/-1sig) points above and below the mean
c     
               if (incr .eq. 1) then
                  plo = pfr
                  phi = i-1
               else
                  plo = i+1
                  phi = pfr
               endif
               do p = plo,phi
                  if (wtau(idx(p)) .ge. tdcf(nbins)) then 
                     j = p
                     goto 50
                 endif
               enddo
50             continue
               p = max0(j-nint(float(j-plo)*0.3413*2.0),plo)
               sigtm(nbins) = tdcf(nbins)-wtau(idx(p))
               p = min0(j+nint(float(phi-j)*0.3413*2.0),phi)
               sigtp(nbins) = wtau(idx(p))-tdcf(nbins)
               pfr = i
               if (pfr .ne. pmax) goto 20
            endif
c If no more points - get out of loop (label 40)
            goto 40
         endif
c     
c Adding another point to the bin...
c     
         if (wuseda(waidx(p)) .or. wusedb(wbidx(p))) goto 30
         inb = inb+1
         if (inb.gt. minbin) then
            print *,'Not enough work area - increase minbin in MAIN!'
            stop
         endif
         wuseda(waidx(p)) = .true.
         wusedb(wbidx(p)) = .true.
         tij = wtau(p)
         tdcf(nbins) = tdcf(nbins)+tij
         vtau(inb,nbins) = tij
         wai(inb,nbins) = waidx(p)
         wbj(inb,nbins) = wbidx(p)
30       continue
      enddo
c     
c Binning is finished
c     
40    continue
      if (.not. (autocf .or. unsmpl .or. incr .eq. 1)) then
         pfr = np/2+1
         pmax = np
         incr = 1
         nnegtv = nbins
         goto 20
      endif
c     
c If CCF (and NOT uniform sampling): Sort the bins into increasing
c chronological order: The nnegtv negative bins are at the beginning but at 
c reverse order.
c     
      if (.not. (autocf .or. unsmpl)) then
         do i = 1,nnegtv/2
            j = nnegtv+1-i
            inb = inbin(i) 
            inbin(i) = inbin(j)
            inbin(j) = inb
            tij = tdcf(i)
            tdcf(i) = tdcf(j)
            tdcf(j) = tij
            tij = sigtp(i)
            sigtp(i) = sigtp(j)
            sigtp(j) = tij
            tij = sigtm(i)
            sigtm(i) = sigtm(j)
            sigtm(j) = tij
            do p = 1,max0(inbin(i),inbin(j))
               tij = vtau(p,i)
               vtau(p,i) = vtau(p,j)
               vtau(p,j) = tij
               inb = wai(p,i)
               wai(p,i) = wai(p,j)
               wai(p,j) = inb
               inb = wbj(p,i)
               wbj(p,i) = wbj(p,j)
               wbj(p,j) = inb
            enddo
         enddo
      endif
c     
c Freeing the allocated work areas
c     
c-malloc(
c      call free(iwaidx)
c      call free(iwbidx)
c      call free(iwtau)
c      call free(iidx)
c-malloc)
c     
c Emergency exit - something is very wrong...
c     
      if (nbins .eq. 0) then
         print *,'alcbin: No bin contains enough points - stopping!'
         stop
      endif
      return
      end
c/////////////////////////////////////////////////////////////////////////////
c Calculating the discrete correlation function.
c POSITIVE lag values mean b lags after a.
c This implementation requires rather big work areas in the interest of a 
c faster algorithm (and is therefore suitable for Monte Carlo simulations). 
c     
      subroutine clcdcf(ta,a,erra,na,
c                       i  i i    i
     >                  tb,b,errb,nb,
c                       i  i i    i
     >                  autocf,no0lag,unsmpl,used,unused,nwdata,
c                       i      i      i      o    o      i/o    
     >                  minpts,tdcf,sigtm,sigtp,dcf,sigdcm,sigdcp,
c                       i/o    o    o     o     o   o      o      
     >                  inbin,minbin,nbins,mbins,
c                       o     i      i/o   i     
     >                  wta,wtb,wa,wb,werra,werrb,wai,wbj,wuseda,wusedb,
c                       w   w   w  w  w     w     w   w   w      w
     >                  vtau)
c                       w
      implicit none
      save
      integer i,j,k,n,used
      integer ibin
      integer na,nb,nbins,minbin,mbins,unused,minpts
      real ta(na),tb(nb)
      real a(na),b(nb)
      real erra(na),errb(nb)
      real xa,xb
      real wa(minbin),wb(minbin)
      real wta(minbin),wtb(minbin)
      real werra(minbin),werrb(minbin)
      integer wai(minbin,mbins),wbj(minbin,mbins)
      real vtau(minbin,mbins)
      real tdcf(mbins),sigtm(mbins),sigtp(mbins),dcf(mbins)
      real sigdcm(mbins),sigdcp(mbins)
      integer inbin(mbins)
      real expa,expb,vara,varb,vnorm,expbin,varbin,z,expz,sigz
      real sqrt,tanh,fishe,fishs,eps
      logical no0lag,nwdata,autocf,unsmpl
      logical wuseda(na),wusedb(nb)
      parameter (eps=1e-7)
c     
c If new data (i.e. NOT another Monte Carlo run) - allocate pairs to bins
c     
      if (nwdata) then
c Allocating the lags to the bins
         call alcbin(ta,na,tb,nb,wuseda,wusedb,no0lag,autocf,minpts,
     >               tdcf,sigtm,sigtp,wai,wbj,inbin,nbins,mbins,minbin,
     >               unsmpl,
     >               vtau)
c Counting the unused points
         used = 0
         do ibin = 1,nbins
            used = used + inbin(ibin)
         enddo
         if (autocf) then
            if (no0lag) then 
               unused = na*(na-1)/2-used
            else
               unused = na*na/2-used
            endif
         else
            unused = na*nb-used
         endif
      endif
c     
c After allocating pairs to bins: calculating the dcf
c     
      do ibin = 1,nbins
c Collecting the points of the bin
         n = inbin(ibin)
         do k = 1,n
            i = wai(k,ibin)
            j = wbj(k,ibin)
            wta(k) = ta(i)
            wtb(k) = tb(j)
            wa(k) = a(i)
            wb(k) = b(j)
            werra(k) = erra(i)
            werrb(k) = errb(j)
         enddo
         expa = 0.0
         expb = 0.0
         vara = 0.0
         varb = 0.0
         do i = 1,n
            xa = wa(i)
            xb = wb(i)
            expa = expa+xa
            expb = expb+xb
            vara = vara+xa**2
            varb = varb+xb**2
         enddo
         expa = expa/n
         expb = expb/n
         vara = (vara-n*expa**2)/(n-1.)
         varb = (varb-n*expb**2)/(n-1.)
c     
         expbin = 0.0
         varbin = 0.0
c If normalization factor is 0 ...
         vnorm = vara*varb
         if (vnorm .le. 0.0) then
            vnorm = 0.0
            expbin = 0.0
         else
            vnorm = sqrt(vnorm)
            do i = 1,n
               expbin = expbin + (wa(i)-expa)*(wb(i)-expb)
            enddo
c Dividing by (n-1) for an unbiased estimator of the correlation coefficient
c cf Barlow / Statistics, p. 80
            expbin = expbin/vnorm/(n-1.)
         endif
         dcf(ibin) = expbin
c     
c Calculating the +/- 1 Sigma limits from Fisher's z
c     
         if      (expbin .ge.  1.0-eps) then
            expbin =  1.0-eps
         else if (expbin .le. -1.0+eps) then 
            expbin = -1.0+eps
         endif
c     
c NOTE: This error estimation is by "bootstrapping": fishe & fishs give
c       the true E(z) and S(z) when the TRUE correlation coefficient is 
c       given. We are using the empirical r itself, similarily to the
c       common poissonian estimate of n +/- sqrt(n)
c     
        z = log((1.+expbin)/(1.-expbin))/2.0
        sigz = fishs(expbin,n)
        expz = fishe(expbin,n)
        sigdcm(ibin) = expbin-tanh(expz-sigz)
        sigdcp(ibin) = tanh(expz+sigz)-expbin
c     
      enddo
      if (nwdata) then
         nwdata = .false.
      endif
      return
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c Fisher's small sample approximation for E(z) (Kendall + Stuart Vol. 1 p.391)
c     
      real function fishe(r,n)
      implicit none
      integer n
      real r
c     
      fishe = log((1.+r)/(1.-r))/2.+
     >        r/2./(n-1.)*(1.+(5.+r**2)/4./(n-1.)+
     >                     (11.+2*r**2+3*r**4)/8./(n-1.)**2)
      return
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c Fisher's small sample approximation for s(z) (Kendall + Stuart Vol. 1 p.391)
c     
      real function fishs(r,n)
      implicit none
      integer n
      real r,sqrt
c     
      fishs = 1./(n-1.)*(1.+(4.-r**2)/2./(n-1.)+
     >                   (22.-6*r**2-3*r**4)/6./(n-1.)**2)
      fishs = sqrt(fishs)
      return
      end
c///////////////////////////////////////////////////////////////////////
c     Generating a standard Gaussian distr. using the Box-Muller method
c     The existence of a subroutine rand(x) is assumed.
      real function rndnrm()
      implicit none
      save y1,y2
      real x1,x2,y1,y2,a
      real rand
      real pi2
      parameter (pi2 = 6.283185307)
      logical odd
      data odd/.true./
c
      if (odd) then
         x1 = rand(0)
         x2 = rand(0)
         x2 = x2*pi2
         a =  sqrt(-2.*log(x1))
         y1 = a*cos(x2)
         y2 = a*sin(x2)
         rndnrm = y1
      else 
         rndnrm = y2
      endif 
      odd = .not. odd
      return 
      end
c///////////////////////////////////////////////////////////////////////
      SUBROUTINE srtind (X, N, IPERM, KFLAG)
C***LIBRARY   SLATEC
C***CATEGORY  N6A1B, N6A2B
C***TYPE      SINGLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
C***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
C***AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SRTIND returns the permutation vector IPERM generated by sorting
C   the array X and, optionally, rearranges the values in X.  X may
C   be sorted in increasing or decreasing order.  A slightly modified
C   quicksort algorithm is used.
C
C   IPERM is such that X(IPERM(I)) is the Ith value in the rearrangement
C   of X.  IPERM may be applied to another array by calling IPPERM,
C   SPPERM, DPPERM or HPPERM.
C
C   The main difference between SRTIND and its active sorting equivalent
C   SSORT is that the data are referenced indirectly rather than
C   directly.  Therefore, SRTIND should require approximately twice as
C   long to execute as SSORT.  However, SRTIND is more general.
C
C   Description of Parameters
C      X - input/output -- real array of values to be sorted.
C          If ABS(KFLAG) = 2, then the values in X will be
C          rearranged on output; otherwise, they are unchanged.
C      N - input -- number of values in array X to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the value in the original order of the
C              X array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting X in increasing order and sort X also.
C            =  1  means return the permutation vector resulting from
C                  sorting X in increasing order and do not sort X.
C            = -1  means return the permutation vector resulting from
C                  sorting X in decreasing order and do not sort X.
C            = -2  means return the permutation vector resulting from
C                  sorting X in decreasing order and sort X also.
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL X(*)
      INTEGER IPERM(*)
C     .. Local Scalars ..
      REAL R, TEMP
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         print *,'SRTIND: N<1'
         stop 
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         print *,'SRTIND: Illegal sort control parameter'
         stop 
      ENDIF
C
C     Initialize permutation vector
C
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE
C
C     Return if only one value is to be sorted
C
      IF (NN .EQ. 1) RETURN
C
C     Alter array X to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 20 I=1,NN
            X(I) = -X(I)
   20    CONTINUE
      ENDIF
C
C     Sort X only
C
      M = 1
      I = 1
      J = NN
      R = .375E0
C
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   40 K = I
C
C     Select a central element of the array and save it in location L
C
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
C
C     If first element of array is greater than LM, interchange with LM
C
      IF (X(IPERM(I)) .GT. X(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than LM, interchange with LM
C
      IF (X(IPERM(J)) .LT. X(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
C
C        If first element of array is greater than LM, interchange
C        with LM
C
         IF (X(IPERM(I)) .GT. X(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
C
C     Find an element in the second half of the array which is smaller
C     than LM
C
   60 L = L-1
      IF (X(IPERM(L)) .GT. X(LM)) GO TO 60
C
C     Find an element in the first half of the array which is greater
C     than LM
C
   70 K = K+1
      IF (X(IPERM(K)) .LT. X(LM)) GO TO 70
C
C     Interchange these elements
C
      IF (K .LE. L) GO TO 50
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 90
C
C     Begin again on another portion of the unsorted array
C
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
C
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
C
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (X(IPERM(I)) .LE. X(LM)) GO TO 100
      K = I
C
  110 IPERM(K+1) = IPERM(K)
      K = K-1
C
      IF (X(LM) .LT. X(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
C
C     Clean up
C
  120 IF (KFLAG .LE. -1) THEN
         DO 130 I=1,NN
            X(I) = -X(I)
  130    CONTINUE
      ENDIF
C
C     Rearrange the values of X if desired
C
      IF (KK .EQ. 2) THEN
C
C        Use the IPERM vector as a flag.
C        If IPERM(I) < 0, then the I-th value is in correct location
C
         DO 150 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0) THEN
               INDX = ISTRT
               INDX0 = INDX
               TEMP = X(ISTRT)
  140          IF (IPERM(INDX) .GT. 0) THEN
                  X(INDX) = X(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               X(INDX0) = TEMP
            ENDIF
  150    CONTINUE
C
C        Revert the signs of the IPERM values
C
         DO 160 I=1,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE
C
      ENDIF
C
      RETURN
      END
c//////////////////////////////////////////////////////////////////////////////
