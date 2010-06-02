cxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxxx    Integrator for CMB anisotropy, CMB polarization and transfer functions
cxxx    Developed and maintained by Uros Seljak (useljak@princeton.edu) and
cxxx    Matias Zaldarriaga (mzaldarriaga@cfa.harvard.edu)
cxxx
cxxx    Permission to use, copy, distribute and modify without fee or royalty at
cxxx    any tier, this software and its documentation, for any purpose
cxxx    and without fee or royalty is hereby granted, provided that
cxxx
cxxx    1)  you agree to comply with the following copyright
cxxx    notice and statements,
cxxx    including the disclaimer, and that the same appear on ALL copies of the
cxxx    software and documentation, including modifications that you make for
cxxx    internal use or for distribution:
cxxx
cxxx    Copyright 1996, 1997, 1998 (versions up to 3.0) by
cxxx    Harvard-Smithsonian Center for
cxxx    Astrophysics and the Massachusetts Institute of Technology.
cxxx    All rights reserved.
cxxx
cxxx    Copyright 1999, 2000 (version 3.0-4.0) by Princeton University
cxxx    and Institute for Advanced Studies, Princeton. All rights
cxxx    reserved.
cxxx
cxxx    Copyright 2003 (version 4.5) by Princeton University
cxxx    and Harvard University. All rights
cxxx    reserved.
cxxx
cxxx    THIS SOFTWARE IS PROVIDED "AS IS", AND INSTITUTIONS ABOVE MAKE NO
cxxx    REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
cxxx    By way of example, but not limitation,
cxxx    THEY MAKE NO REPRESENTATIONS OR WARRANTIES OF
cxxx    MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
cxxx    THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
cxxx    ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
cxxx
cxxx    2)   Applications are limited to non-commercial and not-for-profit
cxxx    purposes.
cxxx
cxxx    3)   An appropriate acknowledgement is included in any publications
cxxx    based on work conducted with usage of the CMBFAST package or any
cxxx    of its derivatives.
cxxx
cxxx    Boltzmann integrator, initial conditions and recombination routines are
cxxx    based on the COSMICS package of E. Bertschinger. See the LICENSE
cxxx    file of the COSMICS distribution for restrictions on the modification
cxxx    and distribution of this software.
cxxx
cxxx    Recombination routine RECFAST written by D. Scott (UBC), based on
cxxx    calculations in the paper by Seager, Sasselov & Scott
cxxx    (ApJ, 523, L1, 1999).
cxxx    See the license in recfast.f for further details regarding restrictions
cxxx    on the modification and distribution of this software.
cxxx
cxxx    Configure script, Makefile, fitsout.f, params.f, FITS and IARGC options
cxxx    written by Matthias Bartelmann.
cxxx
cxxx    parallel version developed with David Spergel and Alexey Makarov
cxxx
cxxx    dark energy/quintessence routines developed with Rob Caldwell
cxxx
cxxxcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                VERSION 4.5                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        program driver
        implicit double precision(a-h,o-z)
c This is a driving routine that illustrates the use of the program.

        include 'cmbfast.inc'
#ifdef MPIBUILD
        include 'mpif.h'
#endif
        include 'cmbpar.inc'

c       const=7*pi**4/120.
        parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
        double precision ztf(ntfmax)
        character*500 filename
        character*500 filename_eos
        character*500 filejlens
        character*500 parafile
        common /filej/ filejlens
        common /lens/ lensflag
        common /trlens/ntlens

c Output arrays temperature: clts, cltt ; e spectra:
c cles, clet ; b perturbation(only tensor contrubute):
c clbt ; cross correlation: clcs, clct.

        double precision clts(l0max,nnmax),cltt(l0max,nnmax)
        double precision cles(l0max,nnmax),clet(l0max,nnmax)
        double precision clbs(l0max,nnmax),clbt(l0max,nnmax)
        double precision clcs(l0max,nnmax),clct(l0max,nnmax)
        double precision clkk(l0max,nnmax),cltk(l0max,nnmax)

        integer l(lmax),l0

        common /lvalues1/ l,l0,lmo
        save /lvalues1/

        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     &                       ,grhonr
        common /initcase/ initfl
        common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
        integer rcflag
        common /recfl/rcflag

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn
        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     &      itflag,lmaxt
        common /transfer/akmaxt,ztf,nlnkt,ict,ntf

        integer ndyn,nflag_rho

        common /qtens/ rcrat,irt
        common /qparm/ wdyn,ndyn

        common /qstore/ nflag_rho

c       K split
        common /cutk/ aksplit, kcutflag

c       TIMING
c       Timing variables for testing purposes
c       Should be commented in released versions
c       real etime,actual,timeprev,tarray(2)
c       external etime
        integer myid,numprocs


c       TIMING
c       Entering CMBFLAT
c       actual=etime(tarray)

#ifdef MPIBUILD
        CALL MPI_INIT_(ierr)
        CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )          !MPI
        CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
#else
        myid=0
        numprocs=1
#endif

#ifdef IARGC
        if (iargc() .eq. 0) then
           intflg = .true.
        else
           intflg = .false.
           if (iargc() .ne. 1) then
              print *, 'Usage: cmb [name of parameter file]'
#ifdef MPIBUILD
              CALL MPI_FINALIZE_(ierr)
#endif
              stop 1
           end if
           call getarg(1,parafile)
           call params(parafile)
        endif
#else
        intflg = .true.
#endif

        do i=1,nnmax
          dalphansdlnk(i)=0.0
        enddo

        if (intflg) then
           IF(MYID.eq.0) THEN
                write(*,*)'CMB (0), transfer functions (1) or both (2):'
                write(*,*)'If you want the lensed Cls you will need (2)'
                read(*,*)ict
           ENDIF
#ifdef MPIBUILD
           CALL MPI_BCAST_(ict,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        else
           ict=p_ict
        endif

        if (ict.ne.1) then
          if (intflg) then
           IF(MYID.eq.0) THEN
                write(*,*)'Value of lmax, ketamax (e.g. 1500 3000)'
                write(*,*)'Remember to be consistent with the file'
                write(*,*)'in the flat case.'
                read(*,*)lmoin,akmax0
            ENDIF
#ifdef MPIBUILD
           CALL MPI_BCAST_(lmoin,1,MPI_INTEGER,0,
     &                          MPI_COMM_WORLD,IERR)
           CALL MPI_BCAST_(akmax0,1,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
          else
           lmoin=p_lmx
           akmax0=p_emx
          endif

           call initlval(lmoin,akmax0)
        endif

       if (ict.ne.0) then
c output at z=10 saves 50% of time; for CDM models tf remains
c unchanged between z=0 and z=10; for MDM, lambda, nonflat models
c one should integrate until z=0
      if (intflg) then
        IF(MYID.eq.0) THEN
        write(*,*)'Enter tf kmax (h/Mpc), # of k per log. int. (5,5)'
        read(*,*)akmaxt,nlnkt
        write(*,*)'Enter number and redshifts of the output: (1,0)'
        write(*,*)'If more than one tf is requested the redshifts'
        write(*,*)'have to be specified in decreasing order.'
        read(*,*)ntf,(ztf(i),i=1,ntf)
        ENDIF

#ifdef MPIBUILD
        CALL MPI_BCAST_(nlnkt,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(akmaxt,1,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(ntf,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(ztf,ntf,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
      else
        akmaxt=p_kmx
        nlnkt=p_nkl
        ntf=p_ntf
        do itf=1,ntf
           ztf(itf)=p_ztf(itf)
        enddo
      endif

       do itf=1,ntf
        if (intflg) then
           IF(MYID.eq.0) THEN
#ifndef LYAOUTPUTS
                write(*,*)'Enter',itf,'. output
     1   transfer function file'
                read(*,'(a500)') filename
#else
                filename='dev/null' ! it's just a trash bin
                if((abs(ztf(itf)-3.0d0).le.0.01d0).or.
     $               (abs(ztf(itf)-0.0d0).le.0.01d0)) then
                   write(*,*)'Enter',itf,
     $                  '. output transfer function file'
                   read(*,'(a500)') filename
                endif
#endif
         ENDIF
        else
            filename=p_ftf(itf)
        endif
        IF(MYID.eq.0) open(unit=12+itf,file=filename,
     1  form='formatted',status='unknown', err=199)
 199    enddo
#ifdef LYAOUTPUTS
        if(MYID.eq.0) then
           write(*,*)'Enter lyaoutputs file'
           read(*,'(a500)') filename
           open(unit=98,file=filename,form='formatted',
     1  status='unknown')
        endif
#endif

        if (intflg.and.MYID.eq.0) then
        write(*,*)'CDM, baryons, photons, massless neutrinos'
        write(*,*)'and massive neutrinos when applicable'
        endif
       endif
c  Read initial parameters.

CC1 If you want to make loops over models remove the
CC1 following questions ad replace them by the loop giving
CC1 values to all the variables,
CC1 omegab,omegac,omegav,omegan, h0,tcmb,yhe,annunr,
CC1 riflag,optdlss or zri,rif,itflag,nn,an,ant (if tensors
CC1 are wanted),initfl.
CC1 Remember that for nonflat models massive neutrinos, gravity
CC1 waves or initial conditions other than isentropic are not
CC1 available.


       if (intflg) then
        nflag_rho = 0
        IF(MYID.eq.0) THEN
          write(*,*)'constant equation of state with dark energy' 
          write(*,*)'perturbations (1), without (3), table with a,w' 
          write(*,*)'with perturbations (2) or without (4)'
          read(*,*)ndyn
        ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(ndyn,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif

          if(ndyn.eq.1.or.ndyn.eq.3) then
                IF(MYID.eq.0) THEN
                        write (*,*) 'Enter wdyn (e.g -1)'
                        read (*,*) wdyn
                ENDIF

#ifdef MPIBUILD
                CALL MPI_BCAST_(wdyn,1,MPI_DOUBLE_PRECISION,0,
     >          MPI_COMM_WORLD,IERR)
#endif
         if(wdyn.eq.-1.0) ndyn=0

       else if(ndyn.eq.2.or.ndyn.eq.4) then
#ifdef MPIBUILD
                if(myid.lt.10) then
                        write(filename_eos,'(a37,i1)')
     >                  'eos.dat.',myid
                else
                        write(filename_eos,'(a37,i2)')
     >                          'eos.dat.',myid
                        endif
#else
        write(*,*)'Enter filename with a,w (e.g.eos.dat)'
        write(*,*)'needs uniform sampling in a'
        write(*,*)' 0<amin<1e-12, amax>1, max elements nq0=',nq0
        read(*,'(a500)')filename_eos
#endif

                open(unit=99,file=filename_eos,
     <          form='formatted',status='unknown')
          endif
         else
c currently only constant wdyn set up for noninteractive mode
          ndyn=p_quin
          wdyn=p_wdyn
          if(wdyn.eq.-1.0) ndyn=0
          if ((ndyn.eq.0).and.(wdyn.ne.-1.0)) then
             write(*,*)'Inconsistent choice  .... '
             write(*,*)'quin=0 means cosmological constant'
             write(*,*)'You then requested w not equal to -1' 
             write(*,*)'Try again'
             stop
          end if
      endif


#ifdef DIM
        if (intflg) then
         IF(MYID.eq.0) then
         write(*,*)'Enter (1) if you want to include
        1        the 5th dimension'
                write(*,*)' Enter(0) if you want the
        2     usual FRW space-time'
                read(*,*)dimflag
         ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(dimflag,1,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif

        else
         dimflag=p_dim
        endif
#endif

        if (intflg) then
           IF(MYID.eq.0) then
                write(*,*)
     2          'Enter Omega_b, Omega_c, Omega_de,',
     1          'Omega_nu (.05 .3 0.65 0)'
                read(*,*) omegab,omegac,omegav,omegan
            ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(omegab,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(omegac,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(omegav,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(omegan,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
#endif
        else
           omegab=p_par(1)
           omegac=p_par(2)
           omegav=p_par(3)
           omegan=p_par(4)
        endif
        omega=omegab+omegac+omegav+omegan

#ifdef DIM
        IF(MYID.eq.0) THEN
                write(*,*)'Enter (1) if you want to
        1       include the 5th dimension'
                write(*,*)' Enter(0) if you want the
        2            usual FRW space-time'
                read(*,*)dimflag
        ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(dimflag,1,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
        if (dimflag.eq.1) omega=(sqrt(omegav)+sqrt(omegav+omegac
     > +omegab+omegan))**2

#endif

        omegak=1.0d0-omega

2       if (intflg) then
          IF(MYID.eq.0) THEN
           write(*,*)
     2    'Enter H0, Tcmb, Y_He, N_nu(massless), N_nu, g*(massive)',
     3    '(e.g. 65 2.726 0.24 3.04 0 10.75)'
                read(*,*) h0,tcmb,yhe,annur,annunr,gsnunr
          ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(h0,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(tcmb,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(yhe,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(annur,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST_(annunr,1,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
#endif
        else
           h0=p_par(5)
           tcmb=p_par(6)
           yhe=p_par(7)
           annur=p_par(8)
           annunr=p_par(9)
           gsnunr=p_par(10)
        endif

        akmaxt=akmaxt*(h0/100.0d0)
        if (h0.lt.25.d0.or.h0.gt.100.d0.and.MYID.eq.0) then
          write(*,*)
     2      '  Warning: H0 has units of km/s/Mpc. Your value is weird.'
        end if
        if (tcmb.lt.2.7d0.or.tcmb.gt.2.8d0.and.MYID.eq.0) then
          write(*,*)
     2      '  Warning: Tcmb has units of K. Your value is weird.'
        end if

        if (yhe.lt.0.2d0.or.yhe.gt.0.3d0.and.MYID.eq.0) then
          write(*,*)
     2      '  Warning: Y_He is the Helium fraction of baryons.',
     3      '  Your value is weird.'
        end if
        if (annunr.lt.0.or.annunr.gt.3.1.and.MYID.eq.0) then
          write(*,*)
     2      'Warning: N_nu(massive) is strange'
#ifdef MPIBUILD
          CALL MPI_FINALIZE_(ierr)
#endif
        end if
        if (annur.lt.0.or.annur.gt.3.1.and.MYID.eq.0) then
          write(*,*)
     2      'Warning: N_nu(massless) is strange'
#ifdef MPIBUILD
          CALL MPI_FINALIZE_(ierr)
#endif
        end if
        if (annunr.lt.1.and.omegan.gt.0.0.and.MYID.eq.0) then
          write(*,*)
     2      'Warning: N_nu(massive) should be 1, 2, or 3',
     3      'For non zero omegan'
#ifdef MPIBUILD
          CALL MPI_FINALIZE_(ierr)
#endif
        end if

        if (intflg) then
        IF(MYID.eq.0) THEN
                write(*,*)  'Enter 0 for Peebles recombination'
                write(*,*)  'or 1 for recfast'
                read(*,*) rcflag
        ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(rcflag,1,MPI_INTEGER,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
        else
           rcflag=p_rcf
        endif

        if (rcflag.eq.1) then
           IF(MYID.eq.0)write(*,*) 'Calling recfast'
           call recfast(omegab,omegac+omegan,omegav,h0,tcmb,yhe,annur)
           IF(MYID.eq.0) write(*,*)'recfast succesfully completed'
        endif

c  Evaluate general cosmology constants.

c  grho gives the contribution to the expansion rate from: (g) photons,
c  (r) one flavor of relativistic neutrino (2 degrees of freedom),
c  (m) nonrelativistic matter (for Omega=1).  grho is actually
c  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
c  a=tau(Mpc)*adotrad, with a=1 today
c  (Used only to set the initial conformal time.)
        grhom=3.3379d-11*h0*h0
        grhog=1.4952d-13*tcmb**4
        grhor=3.3957d-14*tcmb**4
        if (annunr.ne.0) then
         grhonr=grhog*0.875d0*exp(log(43.d0/11.d0/gsnunr)*4.d0/3.d0)
c massive neutrino mass in eV
         xmnu=omegan*h0*h0/10000.0*(gsnunr/10.75)*93.5/annunr
         IF(MYID.eq.0) write(*,*)'neutrino mass=',xmnu,' eV'
        endif

c  adotrad gives the relation a(tau) in the radiation era:
        adotrad=sqrt((grhog+grhor*annur)/3.0)

        if(ndyn.eq.2.or.ndyn.eq.4) call readtable

        if (intflg) then
        IF(MYID.EQ.0) THEN
        write (*,*) 'Enter 0 for no reionization'
        write (*,*) 'Enter 1 for specified optical depth to lss(xe=1)'
        write (*,*) 'Enter 2 for specified redshift and xe'
        read (*,*) riflag
        ENDIF
#ifdef MPIBUILD
        CALL MPI_BCAST_(riflag,1,MPI_DOUBLE_PRECISION,0,
     >          MPI_COMM_WORLD,IERR)
#endif
        else
           riflag=p_rif
        endif

        zri=0.0d0
        rif=0.0d0
        optdlss=0.0d0
        if (riflag.eq.1) then
           if (intflg) then
           IF(MYID.EQ.0) THEN
           write (*,*) 'Enter optical depth to lss'
           read (*,*) optdlss
           ENDIF
#ifdef MPIBUILD
           CALL MPI_BCAST_(optdlss,1,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
           else
              optdlss=p_opd
           endif
           rif=1.0d0
        end if

        if (riflag.eq.2) then
           if (intflg) then
           IF(MYID.EQ.0) THEN
           write (*,*)'Enter redshift, ionization
     1             fraction(e.g. 50 0.2)'
           read (*,*) zri,rif
           ENDIF
#ifdef MPIBUILD
           CALL MPI_BCAST_(zri,1,MPI_DOUBLE_PRECISION,0,
     >                                  MPI_COMM_WORLD,IERR)
           CALL MPI_BCAST_(rif,1,MPI_DOUBLE_PRECISION,0,
     >                                  MPI_COMM_WORLD,IERR)
#endif
           else
              zri=p_red
              rif=p_frc
           endif
           zristp=0.07d0*zri-1.0d0
           if (zristp.lt.0.0) zristp=0.0d0
        end if
#ifdef MPIBUILD
        CALL MPI_BARRIER_(MPI_COMM_WORLD,IERR)
#endif
        if (ict.ne.1) then
           if (intflg) then
             IF(MYID.eq.0) THEN
              write (*,*) 'Enter 0 for scalar alone,'
              write(*,*) '1 for tensor+scalar 2 for tensors alone'
              write(*,*)'3 for scalar (k<k*) and 4 for scalar (k>k*)'
              read(*,*)itflag
             ENDIF
#ifdef MPIBUILD
             CALL MPI_BCAST_(itflag,1,MPI_INTEGER,0,
     1            MPI_COMM_WORLD,IERR)
#endif
           else
              itflag=p_tfl
           endif

c       K splitting
           aksplit=1.0d0
           kcutflag=0
           if (itflag.eq.3) then
              itflag=0
              kcutflag=1
           end if
           if (itflag.eq.4) then
              itflag=0
              kcutflag=-1
           end if

           if (itflag.ne.2) then
              if (intflg) then
              IF(MYID.EQ.0) THEN
              write(*,*)'number and values of scal. spectral index n,'
              write(*,*)'and its running alpha_n (1,1,0)'
              read(*,*)nn,(an(i),i=1,nn),(alphans(i),i=1,nn)
c                 ,(dalphansdlnk(i),i=1,nn)
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(nn,1,MPI_INTEGER,0,
     1             MPI_COMM_WORLD,IERR)
              CALL MPI_BCAST_(an,nn,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
              CALL MPI_BCAST_(alphans,nn,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
c              CALL MPI_BCAST_(dalphansdlnk,nn,MPI_DOUBLE_PRECISION,0,
c     >                  MPI_COMM_WORLD,IERR)
#endif
              else
                 nn=p_nsm
                 do i=1,nn
                    an(i)=p_scl(i)
                    alphans(i)=p_asc(i)
                 enddo
              endif
           end if
           if (itflag.eq.1) then
              if (intflg) then
              IF(MYID.EQ.0) THEN
              write(*,*)
     2 'Tensor spectral index given by nt=ns-1 (0) or different (1)'
              read(*,*)itn
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(itn,1,MPI_INTEGER,0,
     1             MPI_COMM_WORLD,IERR)
#endif
              else
                 itn=p_itn
              endif

              if (itn.eq.0) then
                 do 123 in=1,nn
                    ant(in)=an(in)-1
                    alphant(in)=alphans(in)
 123             continue
                 if (.not.intflg) then
                    do in=1,nn
                       p_ten(in)=ant(in)
                       p_ate(in)=alphant(in)
                    enddo
                 endif

              else
                 if (intflg) then
                 IF(MYID.EQ.0) THEN
                 write(*,*)'values of tensor spectral indexes'
                 write(*,*)'and its running:'
                 read(*,*)(ant(in),in=1,nn),(alphant(i),i=1,nn)
                 ENDIF
#ifdef MPIBUILD
                 CALL MPI_BCAST_(ant,nn,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
#endif
                 else
                    do in=1,nn
                       ant(in)=p_ten(in)
                       alphant(in)=p_ate(in)
                    enddo
                 endif
              endif
           endif
           if (itflag.eq.2) then
              if (intflg) then
              IF(MYID.EQ.0) THEN
              write(*,*)'number, values and running of
     1        tens. spec. ind. (1,0,0)'
              read(*,*)nn,(ant(in),in=1,nn),(alphant(i),i=1,nn)
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(nn,1,MPI_INTEGER,0,
     1             MPI_COMM_WORLD,IERR)
              CALL MPI_BCAST_(ant,nn,MPI_DOUBLE_PRECISION,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
           else
              nn=p_nsm
              do in=1,nn
                 ant(in)=p_ten(in)
                 alphant(in)=p_ate(in)
              enddo
           endif
              do in=1,nn
                 rat(in)=1.0d0
              enddo
           end if
           if (itflag.eq.1) then
              if (intflg) then
              IF(MYID.EQ.0) THEN
              write(*,*)'ratio of tensor to
     1                scalar quadrupole given by'
              write(*,*)'7(1-n_S) (0) or different (1):'
              read(*,*)irt
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(irt,1,MPI_INTEGER,0,
     1             MPI_COMM_WORLD,IERR)
#endif
              else
                 irt=p_irt
              endif

              if (irt.eq.0) then
                 do in=1,nn
                    rat(in)=7.0d0*(1.0d0-an(in))
                 enddo
                 if (.not.intflg) then
                    do in=1,nn
                       p_rat(in)=rat(in)
                    enddo
                 endif
              else
                 if (intflg) then
                 IF(MYID.EQ.0) THEN
                 write(*,*)'values of T/S:'
                 read(*,*)(rat(in),in=1,nn)
                 ENDIF
#ifdef MPIBUILD
                 CALL MPI_BCAST_(rat,nn,MPI_DOUBLE_PRECISION,0,
     >                          MPI_COMM_WORLD,IERR)
#endif

              else
                 do in=1,nn
                    rat(in)=p_rat(in)
                 enddo
              endif
              endif
           end if
           if (itflag.eq.0) then
              do in=1,nn
                 rat(in)=1.0d0
              enddo
           endif

           lensflag=0
           if ((itflag.ne.2).and.(ict.eq.2)) then
              if (intflg) then
              IF(MYID.EQ.0) THEN
              write(*,*)'Enter (0) unlensed Cls only'
              write(*,*)'Enter (1) lensed Cls, linear evolution'
              write(*,*)'Enter (2) lensed Cls, non-linear evolution'
              if (itflag.eq.1) then
                 write(*,*)'ONLY SCALAR Cls ARE LENSED'
              end if
              read(*,*)lensflag
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(lensflag,1,MPI_INTEGER,0,
     >                  MPI_COMM_WORLD,IERR)
#endif
              else
                 lensflag=p_glf
              endif
           end if
           if ((lensflag.ne.0).and.(ict.eq.0)) then
              IF(MYID.eq.0) THEN
              write(*,*)'You did not request the transfer function'
              write(*,*)'calculation needed to do the lensing'
              write(*,*)'you will have to start again'
              ENDIF
#ifdef MPIBUILD
              CALL MPI_FINALIZE_(ierr)
#endif
              stop
           else
c The last requested transfer function will be used for the lensing
              ntlens=ntf
           end if
           IF(MYID.eq.0) THEN
           if (itflag.ne.2) then
              if (intflg) then
              write(*,*)'Enter output filename for SCALAR cl'
              if (lensflag.ne.0) then
                 write(*,*)'If lensing was requested this will'
                 write(*,*)'be the unlensed power spectrum'
              end if
              read(*,'(a500)')filename
              else
                 filename=p_fcl
              endif

              open(unit=8,file=filename,
     &             status='unknown',form='formatted')
              rewind 8

              if (lensflag.ne.0) then
                 if (intflg) then
                 write(*,*)'Enter output filename for
     1                   LENSED SCALAR cl'
                 read(*,'(a500)')filename
                 else
                    filename=p_fgl
                 endif
                 open(unit=10,file=filename,
     &            status='unknown',form='formatted')
                 rewind 10

                 if (intflg) then
                 write(*,*)'Enter file with bessel functions'
                 write(*,*)' for lensing calculation (jlens.dat)'
                 read(*,'(a500)')filejlens
                 else
                    filejlens=p_jgl
                 endif
              end if
           end if

           if (itflag.ne.0) then
              if (intflg) then
              write(*,*)'Enter output filename for TENSOR cl'
              read(*,'(a500)')filename
              else
                 filename=p_tcl
              endif
              open(unit=9,file=filename,
     &             status='unknown',form='formatted')
              rewind 9
           end if
        ENDIF
        else
           itflag=0
           nn=0
           lensflag=0
        end if

#ifdef FITS
        if (intflg) then
              IF(MYID.eq.0) THEN
                        write (*,*) 'Enter FITS file name'
                        read (*,'(a500)') p_fts
              ENDIF
        end if
#endif

        if (itflag.ne.2) then
           if (abs(omegak).gt.1.0d-3) then
              if (intflg) then
              IF(MYID.eq.0) THEN
                write (*,*) 'Enter initial conditions'
                write (*,*) '1= Isentropic (adiabatic)'
                write (*,*) '2= Isocurvature CDM'
                write (*,*) '3= Isocurvature baryon'
                write (*,*) '4= Isocurvature seed conditions'
                read (*,*) initfl
              ENDIF
#ifdef MPIBUILD
              CALL MPI_BCAST_(initfl,1,MPI_INTEGER,0,
     1             MPI_COMM_WORLD,IERR)
#endif

              else
                 initfl=p_inf
              endif
           else
              if (intflg) then
              IF(MYID.eq.0) THEN
                write (*,*) 'Enter initial conditions'
                write (*,*) '1= Isentropic (adiabatic)'
                write (*,*) '2= Isocurvature CDM'
                write (*,*) '3= Isocurvature baryon'
                write (*,*) '4= Isocurvature seed conditions'
                read (*,*) initfl
               ENDIF
#ifdef MPIBUILD
               CALL MPI_BCAST_(initfl,1,MPI_INTEGER,0,
     1              MPI_COMM_WORLD,IERR)
#endif

              else
                 initfl=p_inf
              endif
           end if
        else
           initfl=1
        end if

        if (intflg) then
           p_frc=rif
           p_glf=lensflag
           p_ict=ict
           p_inf=initfl
           p_irt=irt
           p_itn=itn
           p_kmx=akmaxt
           p_nkl=nlnkt
           p_nsm=nn
           p_ntf=ntf
           p_opd=optdlss
           p_par(1)=omegab
           p_par(2)=omegac
           p_par(3)=omegav
           p_par(4)=omegan
           p_par(5)=h0
           p_par(6)=tcmb
           p_par(7)=yhe
           p_par(8)=annur
           p_par(9)=annunr
           p_rat(1)=rat(1)
           p_red=zri
           p_rcf=rcflag
           p_rif=riflag
           p_scl(1)=an(1)
           p_ten(1)=ant(1)
           p_tfl=itflag
           p_ztf(1)=ztf(1)
        end if
        if (.not.intflg) call getpar(.false.)

c     TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(*,*)actual-timeprev,' END INPUT, NOW MAIN'
c       write(50,*)'PARAMETERS'
c       write(50,*)omegab,omegac,omegav,omegan
c       write(50,*)h0,tcmb,yhe,annur,annunr
c       write(50,*)riflag,zri,rif,optdlss

c call main subroutine

CC2 If you are interested in constructing a driver for
CC2 only the flat or nonflat models code, leave the appropiate
CC2 piece of following if statement. Then enter the appropiate
CC2 parameters for that case.

        do in=1,nnmax
           do il=1,l0max
              clts(il,in)=0.0d0
              cltt(il,in)=0.0d0
              cles(il,in)=0.0d0
              clet(il,in)=0.0d0
              clbs(il,in)=0.0d0
              clbt(il,in)=0.0d0
              clcs(il,in)=0.0d0
              clct(il,in)=0.0d0
              clkk(il,in)=0.0d0
              cltk(il,in)=0.0d0
           end do
        end do

       if (abs(omegak).le.1.0d-3) then
CC1 initjl reads the tables of spherical Bessel functions
CC1 so if you are making a loop over models you should
CC1 move this call outside the loop.
          if (ict.ne.1) then

             if (intflg) then
#ifdef MPIBUILD
                IF(MYID.eq.0) THEN
c               write(*,*)'Enter input filename for jl'
c               read(*,'(a500)')filename
                write(*,*) ' reading in jl data from
     1                  Data/jl.dat.myid'
                ENDIF
                if(myid.lt.10) then
                write(filename,'(a37,i1)')
     >                  'Data0/jl.dat.',myid
                else
                write(filename,'(a37,i2)')
     >                  'Data0/jl.dat.',myid
                endif
#else
                write(*,*)'Enter input filename for jl'
                read(*,'(a500)')filename
#endif
             else
                filename=p_fjl
             endif
             call initjl(filename)
          end if
          call cmbflat(clts,cltt,cles,clet,
     1         clbt,clcs,clct,clkk,cltk)
        else
           if (ict.ne.1) then
              if (intflg) then
#ifdef MPIBUILD
                IF(MYID.eq.0) THEN
                write(*,*) ' reading in ujl
     1                  data from Data/ujl.dat.myid'
                ENDIF
                if(myid.lt.10) then
                write(filename,'(a38,i1)')
     >                  'ujl.dat.',myid
                else
                write(filename,'(a38,i2)')
     >                  'ujl.dat.',myid
                endif
#else
                write(*,*)'Enter input filename for ujl'
                read(*,'(a500)')filename
#endif
              else
                 filename=p_ujl
              endif
              call initujl0(filename) !modified for MPI usage
           end if
           call cmbopen(clts,cltt,cles,clet,
     1          clbt,clcs,clct,clkk,cltk)
        end if

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(*,*)actual-timeprev,' END MAIN, NOW COBE'


        IF(MYID.eq.0) THEN
c For many purposes it is useful not to call
c COBE normalization so that there is a direct relation
c between the initial power spectrum and the CMB.

#ifdef UNNORM
c       Not calling COBE normalization
#else
        call COBEnormalize(clts,cltt,cles,clet,clbt,
     &  clcs,clct,clkk,cltk)
#endif

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(*,*)actual-timeprev,' END COBE, NOW OUTPUT'


        call output(clts,cltt,cles,clet,clbt,
     &  clcs,clct,itflag,lmo)

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(*,*)actual-timeprev,' END OUTPUT, NOW LENSING'


        if (lensflag.ne.0) then
           call lensing(clts,cles,clbs,clcs,clkk,cltk)
        end if

        if (itflag.ne.2) close(8)
        if (itflag.ne.0) close(9)

#ifdef FITS
        call fitsout(clts,cles,clbt,clcs,itflag,lmo)
#endif


        ENDIF
c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(*,*)actual-timeprev,' END LENSING, NOW END'

CC1 End your Loop over models here.
#ifdef MPIBUILD
       CALL MPI_FINALIZE_(ierr)
#endif
       stop
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
