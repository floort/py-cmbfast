c --------------------------------------------------------------------
c fitsout:
c write transfer function as FITS table
c M. Bartelmann, Jan. 18, 1999
c E. Hivon, Feb. 4, 2003
c  -edited to make FITS output in the CMBFAST (EB)
c   normalization convention instead of CG (a la Kamionkowsi et al.)
c   since Healpix 1.2 now matches CMBFAST
c  -added FITS keyword POLNORM = CMBFAST
c  -changed version number in FITS file to 4.01
c  -let cftisio write date in Y2K format
c --------------------------------------------------------------------
      subroutine fitsout(clts,cles,clbt,clcs,itflag,lmo)
      include 'cmbfast.inc'
      include 'cmbpar.inc'

      parameter (pi=3.141592654)
c       parameter (nh=50)
      parameter (nh=55)
      double precision clts(l0max,nnmax)
      double precision cles(l0max,nnmax)
      double precision clbt(l0max,nnmax)
      double precision clcs(l0max,nnmax)
      real cl(0:l0max,4)
      character*80 hd(nh)

      lmx=lmo-300
      do i=0,1
         do j=1,4
            cl(i,j)=0.0
         enddo
      enddo

      tcmb=p_par(6)
      tcmb2=tcmb*tcmb
      do i=2,lmx
         f=real(i*(i+1))/2.0/pi

         cl(i,1)=real(clts(i,1))*tcmb2/f
c        CG convention was used by Healpix 1.1
c          cl(i,2)=real(cles(i,1))*tcmb2/f/2.0
         cl(i,2)=real(cles(i,1))*tcmb2/f
         cl(i,3)=0.0
c        CG convention was used by Healpix 1.1
c         cl(i,4)=real(clcs(i,1))*tcmb2/f/sqrt(2.0)
         cl(i,4)=real(clcs(i,1))*tcmb2/f
      end do

      if (itflag.ne.0) then
         do i=2,lmx
            f=real(i*(i+1))/2.0/pi

c           CG convention was used by Healpix 1.1
c            cl(i,3)=real(clbt(i,1))*tcmb2/f/2.0
            cl(i,3)=real(clbt(i,1))*tcmb2/f
         end do
      end if

      c2=cl(2,1)
      qr=sqrt(2.5*c2)*1.0e6

      print '(" Q_rms-ps = ",e12.4," micro K")', qr

      call makeHeader(nh,hd,qr)
      call asctab(cl,lmx,4,hd,nh,p_fts)

      return
      end
c --------------------------------------------------------------------
c     clout    = power spectra with l in (0:lmo)
c     ncl      = number of spectra
c     header   = FITS header to be put on top of the file
c     nlheader = number of lines of the header
c     filename = FITS output file name
c --------------------------------------------------------------------
      subroutine asctab(clout,lmo,ncl,hd,nh,file)
      implicit none
      integer l0max,lmax,nnmax,nk0,nstep0
      include 'cmbfast.inc'
      integer nclmax
      parameter (nclmax=20)

      integer lmo,ncl,nh,i
      real clout(0:l0max,1:ncl)
      character*(*) file
      character*80 hd(1:nh)

      integer status,unit,blocksize,bitpix,nspace
      integer naxis,naxes(1),nrows,rowlen,frow,felem,colnum
      integer tfields,tbcol(nclmax)
      logical simple,extend
      character*16 tform(nclmax),ttype(nclmax),tunit(nclmax)
      character*16 extname
c nor longer used [MR]
c      character*80 old_dat
c      character*10 fulldat
c      character*80 comment

      status=0
      blocksize=1
      bitpix=32
      naxis=0
      naxes(1)=0
      simple=.true.
      extend=.true.

      call ftgiou(unit,status)
      call ftinit(unit,file,blocksize,status)
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      do i=1,3
         if (hd(i).ne.' ') call ftprec(unit,hd(i),status)
      end do
      call ftpdat(unit,status)
      if (status.gt.0) call errorMessage(status)

c     commented out. cftisio is now Y2K compliant [EH]
c       ! update the date (format ccyy-mm-dd) to conform with FITS
c       ! requirements. Will need replacing when cfitsio is compliant
c       call ftgkys(unit,'DATE',old_dat,comment,status)
c       fulldat = '20'//old_dat(7:8)//'-'//old_dat(4:5)//'-'//old_dat(1:2)
c       comment = 'FITS file creation date ccyy-mm-dd'
c       call ftukys(unit,'DATE',fulldat,comment,status)
      if (status.gt.0) call errorMessage(status)

      call ftcrhd(unit,status)
      if (status.gt.0) call errorMessage(status)
      nrows=lmo+1
      tfields=ncl
      ttype(1)='Temperature C_l'
      ttype(2)='E-mode C_l'
      ttype(3)='B-mode C_l'
      ttype(4)='T-E cross corr.'
      do i=1,ncl
         tform(i)='E15.7'
         tunit(i)='Kelvin-squared'
      end do
      extname=' '
      nspace=1
      call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)
      if (status.gt.0) call errorMessage(status)
      call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit,
     &     extname,status)
      if (status.gt.0) call errorMessage(status)
      do i=3,nh
         if (hd(i).ne.' ') call ftprec(unit,hd(i),status)
      end do
      if (status.gt.0) call errorMessage(status)
      frow=1
      felem=1
      do colnum=1,tfields
         call ftpcle(unit,colnum,frow,felem,nrows,clout(0,colnum),
     &        status)
      end do
      if (status.gt.0) call errorMessage(status)
      call ftclos(unit,status)
      if (status.gt.0) call errorMessage(status)
      return
      end

c --------------------------------------------------------------------
c create header
c --------------------------------------------------------------------

      subroutine makeHeader(nh,hd,qr)
      implicit none
      include 'cmbpar.inc'
      integer i,nh
      real qr
      character*80 hd(nh),card

      do i=1,nh
         hd(i)=' '
      end do

      call headerLine('COMMENT --------------------',hd(1))
      call headerLine('CREATOR = CMBFAST_4.4',hd(2))
      call headerLine('COMMENT --------------------',hd(3))
      call headerLine('COMMENT Cosmological parameters',hd(4))
      call headerLine('COMMENT --------------------',hd(5))
      write (card,
     &     '("wdyn = ",f8.2," / Equation of state parameter")')
     &     p_wdyn
      call headerLine(card,hd(6))
      write (card,
     &     '("OMEGAB = ",f8.2," / Omega in baryons")')
     &     p_par(1)
      call headerLine(card,hd(7))
      write (card,
     &     '("OMEGAC = ",f8.2," / Omega in dark matter")')
     &     p_par(2)
      call headerLine(card,hd(8))
      write (card,
     &     '("OMEGAV = ",f8.2," / Omega in cosmological constant")')
     &     p_par(3)
      call headerLine(card,hd(9))
      write (card,
     &     '("OMEGAN = ",f8.2," / Omega in neutrinos")')
     &     p_par(4)
      call headerLine(card,hd(10))
      write (card,
     &     '("HUBBLE = ",f8.2," / Hubble constant in km/s/Mpc")')
     &     p_par(5)
      call headerLine(card,hd(11))
      write (card,
     &     '("NNUNR = ",f8.2," / number of massive neutrinos")')
     &     p_par(9)
      call headerLine(card,hd(12))
      write (card,
     &     '("NNUR = ",f8.2," / number of massless neutrinos")')
     &     p_par(8)
      call headerLine(card,hd(13))
      write (card,
     &     '("TCMB = ",f8.4," / CMB temperature in Kelvin")')
     &     p_par(6)
      call headerLine(card,hd(14))
      write (card,
     &     '("HELFRACT = ",f8.2," / Helium fraction")')
     &     p_par(7)
      call headerLine(card,hd(15))
      write (card,
     &     '("OPTDLSS = ",f8.2," / reionisation optical depth")')
     &     p_opd
      call headerLine(card,hd(16))
      write (card,
     &     '("IONFRACT = ",f8.2," / ionisation fraction")')
     &     p_frc
      call headerLine(card,hd(17))
      write (card,
     &     '("ZREION = ",f8.2," / reionisation redshift")')
     &     p_red
      call headerLine('COMMENT --------------------',hd(18))
      call headerLine('COMMENT Modes',hd(19))
      call headerLine('COMMENT --------------------',hd(20))
      write (card,
     &     '("NMODES = ",i8," / number of perturbation modes")')
     &     p_nsm
      call headerLine(card,hd(21))
      write (card,
     &     '("SCLIDX1 = ",f8.2," / scalar power-spectrum index")')
     &     p_scl(1)
      call headerLine(card,hd(22))
      write (card,
     &     '("TENIDX1 = ",f8.2," / tensor power-spectrum index")')
     &     p_ten(1)
      call headerLine(card,hd(23))
      write (card,
     &     '("TENRAT1 = ",f8.2," / tensor-to-scalar ratio")')
     &     p_rat(1)
      call headerLine(card,hd(24))
      call headerLine('COMMENT --------------------',hd(25))
      call headerLine('COMMENT Transfer-function parameters',hd(26))
      call headerLine('COMMENT --------------------',hd(27))
      write (card,
     &     '("MAXK = ",f8.2," / Maximum k in transfer function")')
     &     p_emx
      call headerLine(card,hd(28))
      write (card,
     &     '("MAXKMPC = ",f8.2," /Mpc")')
     &     p_kmx
      call headerLine(card,hd(29))
      write (card,
     &     '("NKPERLOG = ",i8," / number of k per log interval")')
     &     p_nkl
      call headerLine(card,hd(30))
      write (card,
     &     '("NTRANSF = ",i8," / number of transfer functions")')
     &     p_ntf
      call headerLine(card,hd(31))
      write (card,
     &     '("ZTRANS1 = ",f8.2," / transfer-function redshift")')
     &     p_ztf(1)
      call headerLine(card,hd(32))
      write (card,
     &     '("MAXMULT = ",i8," / maximum l")')
     &     p_lmx
      call headerLine(card,hd(33))
      call headerLine('COMMENT --------------------',hd(34))
      call headerLine('COMMENT Flags',hd(35))
      call headerLine('COMMENT --------------------',hd(36))
      write (card,
     &     '("RECFLAG = ",i8," / recombination flag")')
     &     p_rcf
      call headerLine(card,hd(37))
      write (card,
     &     '("RIONFLAG = ",f8.2," / reionisation flag")')
     &     p_rif
      call headerLine(card,hd(38))
      write (card,
     &     '("CALCTYPE = ",i8," / calculation type")')
     &     p_ict
      call headerLine(card,hd(39))
c      write (card,
c     &     '("PRECFLAG = ",i8," / precision flag")')
c     &     p_prc
c      call headerLine(card,hd(40))
      write (card,
     &     '("INITFLAG = ",i8," / initial-condition flag")')
     &     p_inf
      call headerLine(card,hd(40))
      write (card,
     &     '("TRATFLAG = ",i8," / tensor-to-scalar ratio flag")')
     &     p_irt
      call headerLine(card,hd(41))
      write (card,
     &     '("MODEFLAG = ",i8," / selection of modes")')
     &     p_tfl
      call headerLine(card,hd(42))
      write (card,
     &     '("TIDXFLAG = ",i8," / tensor-mode index flag")')
     &     p_itn
      call headerLine(card,hd(43))
      write (card,
     &     '("LENSFLAG = ",i8," / gravitational-lensing flag")')
     &     p_glf
      call headerLine(card,hd(44))
      call headerLine('COMMENT --------------------',hd(45))
      call headerLine('COMMENT COBE Normalisation',hd(46))
      call headerLine('COMMENT --------------------',hd(47))
      write (card,
     &     '("QRMS-PS = ",e12.4," / Q_rms-PS in micro-Kelvin")')
     &     qr
      call headerLine(card,hd(48))
c     add keyword on polarisation convention
      call headerLine('COMMENT --------------------',hd(49))
      call headerLine('COMMENT Normalisation convention',hd(50))
      call headerLine('COMMENT --------------------',hd(51))
      call headerLine(
     &     'POLNORM = CMBFAST / CMBFAST (EB) convention for '
     &     // 'polar. spectra'
     &     ,   hd(52))
      call headerLine('COMMENT --------------------',hd(53))

      return
      end

c --------------------------------------------------------------------
c return formatted header line
c --------------------------------------------------------------------

      subroutine headerLine(line,hdLine)
      implicit none
      integer hdtype,status
      character*(*) line
      character*80 hdLine

      hdtype=0
      status=0
      call ftgthd(line,hdLine,hdtype,status)
      if (status.gt.0) call errorMessage(status)
      end

c --------------------------------------------------------------------
c print error message and stop
c --------------------------------------------------------------------

      subroutine errorMessage(status)
      implicit none
      integer status
      character*30 errmsg

      call ftgerr(status,errmsg)
      print '("FITS error ",i3,": ",a)', status,errmsg
      stop
      end
