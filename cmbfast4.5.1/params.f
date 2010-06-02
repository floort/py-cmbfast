c --------------------------------------------------------------------
c params:
c parse parameter file and assign parameters in common block [cmbpar]
c M. Bartelmann, Jan. 18, 1999
c --------------------------------------------------------------------

      subroutine params(filename)
      implicit none
      include 'cmbpar.inc'
      integer i,j,k,scan,ind
      character*500 line,name,value,adjustl
      character*(*) filename

      external adjustl,scan

      call setpar()

      open (1,file=filename,status='old')
    1 read (1,'(a)',end=2) line
      i=scan(line,'=')
      if (i.gt.0 .and. line(1:1).ne.'#') then
         name=adjustl(line(1:i-1))
         value=adjustl(line(i+1:))
         j=scan(name,'(')
         if (j.gt.0) then
            k=scan(name,')')
            if (k.eq.0) stop 'error 1 in reader'
            read (name(j+1:k-1),*) ind
            if (ind.ge.10) stop 'error 2 in reader'
            if (name(1:j-1).eq.'z')   read (value,*) p_ztf(ind)
            if (name(1:j-1).eq.'an')  read (value,*) p_scl(ind)
            if (name(1:j-1).eq.'alphans')  read (value,*) p_asc(ind)
c            if (name(1:j-1).eq.'dalphansdlnk')
c     $           read (value,*) p_das(ind)
            if (name(1:j-1).eq.'ant') read (value,*) p_ten(ind)
            if (name(1:j-1).eq.'alphant')  read (value,*) p_ate(ind)
            if (name(1:j-1).eq.'rat') read (value,*) p_rat(ind)
            if (name(1:j-1).eq.'ftf') p_ftf(ind)=value
         else
            if (name.eq.'initfl')   read (value,*) p_inf
            if (name.eq.'ict')      read (value,*) p_ict
            if (name.eq.'lmo')      read (value,*) p_lmx
            if (name.eq.'akmax0')   read (value,*) p_emx
            if (name.eq.'akmaxt')   read (value,*) p_kmx
            if (name.eq.'nlnkt')    read (value,*) p_nkl
            if (name.eq.'ntf')      read (value,*) p_ntf
            if (name.eq.'omegab')   read (value,*) p_par(1)
            if (name.eq.'omegac')   read (value,*) p_par(2)
            if (name.eq.'omegav')   read (value,*) p_par(3)
            if (name.eq.'omegan')   read (value,*) p_par(4)
            if (name.eq.'h0')       read (value,*) p_par(5)
            if (name.eq.'tcmb')     read (value,*) p_par(6)
            if (name.eq.'yhe')      read (value,*) p_par(7)
            if (name.eq.'annur')    read (value,*) p_par(8)
            if (name.eq.'annunr')   read (value,*) p_par(9)
            if (name.eq.'gsnunr')   read (value,*) p_par(10)
            if (name.eq.'5dim')     read (value,*) p_dim
            if (name.eq.'quin')     read (value,*) p_quin
            if (name.eq.'wdyn')     read (value,*) p_wdyn
            if (name.eq.'rcflag')   read (value,*) p_rcf
            if (name.eq.'riflag')   read (value,*) p_rif
            if (name.eq.'optdlss')  read (value,*) p_opd
            if (name.eq.'zri')      read (value,*) p_red
            if (name.eq.'rif')      read (value,*) p_frc
            if (name.eq.'itflag')   read (value,*) p_tfl
            if (name.eq.'nn')       read (value,*) p_nsm
            if (name.eq.'itn')      read (value,*) p_itn
            if (name.eq.'irt')      read (value,*) p_irt
            if (name.eq.'lensflag') read (value,*) p_glf
            if (name.eq.'fcl')      p_fcl=value
            if (name.eq.'fgl')      p_fgl=value
            if (name.eq.'fjl')      p_fjl=value
            if (name.eq.'tcl')      p_tcl=value
            if (name.eq.'ujl')      p_ujl=value
            if (name.eq.'fts')      p_fts=value
            if (name.eq.'jgl')      p_jgl=value
         end if
      end if
      goto 1

    2 if (p_glf.ne.0) p_ict=2

      return
      end

c --------------------------------------------------------------------
c scan(s,c):
c search for character 'c' in string 's'
c --------------------------------------------------------------------
      function scan(s,c)
      implicit none
      integer i,scan
      character*(*) s
      character*1 c,x

      scan=0
      do i=1,79
         x=s(i:i+1)
         if (x.eq.c) scan=i
      end do

      return
      end

c --------------------------------------------------------------------
c adjustl(s):
c remove leading blanks from string 's'
c --------------------------------------------------------------------
      function adjustl(s)
      implicit none
      character*(*) s
      character*500 t,adjustl

      t=s
      if (t.eq.' ') goto 2

    1 if (t(1:1).ne.' ') goto 2
      t=t(2:)
      goto 1

    2 adjustl=t
      return
      end

c --------------------------------------------------------------------
c setpar:
c set parameters to default values
c M. Bartelmann, Jan. 18, 1999
c --------------------------------------------------------------------
      subroutine setpar()
      implicit none
      include 'cmbpar.inc'

      p_ict = 0
      p_lmx = 1500
      p_emx = 3000.0
      p_kmx = 5.0
      p_nkl = 5

      p_ntf = 1
      p_ztf(1) = 0.0
      p_ftf(1) = 'trans.d'

      p_par(1) = 0.05
      p_par(2) = 0.95
      p_par(3) = 0.0
      p_par(4) = 0.0
      p_par(5) = 50.0
      p_par(6) = 2.726
      p_par(7) = 0.24
      p_par(8) = 3.04
      p_par(9) = 0.0
      p_par(10) = 10.75

      p_quin= 0
      p_wdyn= -1
      p_rcf = 0

      p_rif = 0.0
      p_opd = 0.0
      p_red = 50.0
      p_frc = 0.2

      p_tfl = 0
      p_nsm = 1
      p_scl(1) = 1.0
      p_asc(1) = 0.0
      p_itn = 0
      p_ten(1) = 0.0
      p_ate(1) = 0.0
      p_irt = 0
      p_rat(1) = 0.0

      p_glf = 0
      p_fcl = 'cl_unlensed.d'
      p_fgl = 'cl_lensed.d'
      p_jgl = 'jlens.dat'

      p_tcl = 'cl_tensor.d'

      p_inf = 1

      p_fjl = 'jl.dat'
      p_ujl = 'ujl.dat'
      p_fts = 'cl.fits'

      return
      end

c --------------------------------------------------------------------
c getpar:
c print parameter values
c M. Bartelmann, Jan. 19, 1999
c --------------------------------------------------------------------
      subroutine getpar(all)
      implicit none
      integer i
      logical all
      include 'cmbpar.inc'

      print 9, 'CMB Parameters','--------------'
    9 format(/,a14,/,a14)

      print 1, p_ict,p_lmx,p_emx,p_kmx,p_nkl
    1 format(
     &     "ict      = ",i5,/,
     &     "lmo      = ",i5,"   , akmax0   = ",f8.2,/,
     &     "akmaxt   = ",f8.2,", nlnkt    = ",i5)

      print 2, p_ntf
    2 format("ntf      = ",i5)
      if (all) then
         do i=1,p_ntf
            print 3, i,p_ztf(i),i,p_ftf(i)
    3       format("z(",i1,")     = ",f8.2,", file(",i1,")  = ",a)
         enddo
      else
         do i=1,p_ntf
            print 13, i,p_ztf(i)
   13       format("z(",i1,")     = ",f8.2)
         enddo
      endif

      print 23, p_quin,p_wdyn
   23   format(
     &     "ndyn     = ",i8,", wdyn     = ",f8.2)

      print 4, (p_par(i),i=1,10)
    4 format(
     &     "omegab   = ",f8.2,", omegac   = ",f8.2,/,
     &     "omegav   = ",f8.2,", omegan   = ",f8.2,/,
     &     "h0       = ",f8.2,", tcmb     = ",f8.2,/,
     &     "yhe      = ",f8.2,", annunr   =",f8.2,/,
     &     "annunr   = ",f8.2,", gsnunr   = ",f8.2)

      print 5, p_rcf,p_rif,p_opd,p_red,p_frc
    5 format(
     &     "rcflag   = ",i8,/,
     &     "riflag   = ",f8.2,", optdlss  = ",f8.2,/,
     &     "zri      = ",f8.2,", rif      = ",f8.2)

      print 6, p_tfl,p_nsm,p_itn,p_irt
    6 format(
     &     "itflag   = ",i5,"   , nn       = ",i5,/,
     &     "itn      = ",i5,"   , irt      = ",i5)
      do i=1,p_nsm
         print 7, i,p_scl(i),i,p_asc(i),i,p_ten(i),i,p_ate(i),i,p_rat(i)
    7    format(
     &        "an(",i1,")         = ",f8.2,/,
     &        "alphans(",i1,")    = ",f8.2,/,
c     &        "dalphansdlnk(",i1,")    = ",f8.2,/,
     &        "ant(",i1,")        = ",f8.2,/,
     &        "alphant(",i1,")    = ",f8.2,/,
     &        "rat(",i1,")        = ",f8.2)
      enddo

      if (all) then
         print 8, p_glf,p_inf,p_fcl,p_fgl,p_jgl,p_tcl,p_fjl,p_ujl,p_fts
    8    format(
     &        "lensflag = ",i5,"   , initfl   = ",i5,/,/,
     &        "file (unlensed scalar C_l) = ",a,/,
     &        "file (lensed scalar C_l)   = ",a,/,
     &        "file (jl for lensing)      = ",a,/,
     &        "file (tensor C_l)          = ",a,/,
     &        "file (jl)                  = ",a,/,
     &        "file (ujl)                 = ",a,/,
     &        "file (FITS)                = ",a,/)
      else
         print 18, p_glf,p_inf
   18    format(
     &        "lensflag = ",i5,"   , initfl   = ",i5,/)
      endif

      return
      end
