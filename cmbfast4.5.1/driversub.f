C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (matiasz@sns.ias.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine run_cmbfast(
     &          omega_b,omega_c,omega_v,omega_n,rwdyn,
     &          H_0,t_cmb,Y_he,a_nu_r,a_nu_nr,
     &          tau,num_slope,slope_n,slope_al,ntype,l_max,aketamax,
     &          outfile)

        implicit double precision(a-h,o-z)
c This is a driving routine that illustrates the use of the program.

        real omega_b,omega_c,omega_v,omega_n
        real H_0,t_cmb,Y_he,a_nu_r,a_nu_nr
        real aketamax
        real rwdyn
        include 'cmbfast.inc'

c       parameter (nnmax=10)

        real tau, slope_n(20),slope_al(20)
        integer ntype, type ! 0 = scalar, 1 = tensor (does EITHER one OR the other)
                        ! slope_n refers to either scalar or tensor index, depending on type.
        integer l_max   ! l-value up to which C_l should be computed
        integer num_slope ! number of ns to be calculated must be less than 20


c       const=7*pi**4/120.
        parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
        double precision ztf(nnmax)
        character*500 filename, outfile
        character*500 filejlens
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
        common /initcase/ initfl
        common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
        integer rcflag
        common /recfl/rcflag


        common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn

        common /tensor/ ant(nnmax),rat(nnmax),itflag,lmaxt
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

c       TIMING
c       Enterning CMBFLAT
c       actual=etime(tarray)
c       timeinit=actual


c       FOR PROJECT WITH MAX

        if ((ntype.ne.0).and.(ntype.ne.1).and.(ntype.ne.2)
     $          .and.(ntype.ne.3)) then
           write(*,*)'Invalid type'
           stop
        end if

c       KSPLIT
        aksplit=1.0d0
        kcutflag=0
        if (ntype.eq.0) then
           type = 0
           kcutflag=0
        end if
        if (ntype.eq.1) then
           type = 1
           kcutflag=0
        end if
        if (ntype.eq.2) then
           type = 0
           kcutflag=1
        end if
        if (ntype.eq.3) then
           type = 0
           kcutflag=-1
        end if

c       TIMING
c       Enterning CMBFLAT
c       actual=etime(tarray)

c       THIS DRIVER IS JUST FOR COMPUTING CMB
c       write(*,*)'CMB (0), transfer functions (1) or both (2):'
c       write(*,*)'If you want the lensed Cls you will need (2)'
c       read(*,*)ict
        ict=0

c          write(*,*)'Value of lmax, ketamax (e.g. 1500 3000)'
c          write(*,*)'Remember to be consistent with the file'
c          write(*,*)'in the flat case.'
c          read(*,*)lmoin,akmax0
c
c

        if (aketamax.eq.0.0) then
           akmax0=200.0d0+2.0d0*dble(l_max)
           akmax0=min(akmax0,5000.d0)
        else
           akmax0=dble(aketamax)
        end if
        call initlval(l_max,akmax0)



       if (ict.ne.0) then
c output at z=10 saves 50% of time; for CDM models tf remains
c unchanged between z=0 and z=10; for MDM, lambda, open models
c one should integrate until z=0
        write(*,*)'Enter tf kmax (h/Mpc), # of k per log. int. (5,5)'
        read(*,*)akmaxt,nlnkt
        write(*,*)'Enter number and redshifts of the output: (1,0)'
        write(*,*)'If more than one tf is requested the redshifts'
        write(*,*)'have to be specified in decreasing order.'
        read(*,*)ntf,(ztf(i),i=1,ntf)
        do itf=1,ntf
         write(*,*)'Enter',itf,'. output transfer function file'
         read(*,'(a500)') filename
         open(unit=12+itf,file=filename,form='formatted'
     <,status='unknown')
        enddo
        write(*,*)'output file consists of k(h/Mpc) and tf for '
        write(*,*)'CDM, baryons, photons, massless neutrinos'
        write(*,*)'and massive neutrinos when applicable'
       endif

c  Read initial parameters.

CC1 If you want to make loops over models remove the
CC1 following questions ad replace them by the loop giving
CC1 values to all the variables,
CC1 omegab,omegac,omegav,omegan, h0,tcmb,yhe,annunr,
CC1 riflag,optdlss or zri,rif,itflag,nn,an,ant (if tensors
CC1 are wanted),initfl.
CC1 Remember that for open models massive neutrinos, gravity
CC1 waves or initial conditions other than isentropic are not
CC1 available.

        nflag_rho=0
c       constant equation of state
        ndyn=1
c       wdyn, choose cosmologiacal constan
        wdyn=dble(rwdyn)
        if (wdyn.eq.-1.0)ndyn=0


c1      write(*,*)
c     2  'Enter Omega_b, Omega_c, Omega_v, Omega_nu (e.g. .05 .95 0 0)'
c       read(*,*) omegab,omegac,omegav,omegan
        omegab  = dble(omega_b)
        omegac  = dble(omega_c)
        omegav  = dble(omega_v)
        omegan  = dble(omega_n)
        omega=omegab+omegac+omegav+omegan
        omegak=1.0d0-omega

 2      continue
c2      write(*,*)
c     2    'Enter H0, Tcmb, Y_He, N_nu(massless), N_nu(massive)',
c     3    '(e.g. 50 2.726 0.24 3.04 0)'
c       read(*,*) h0,tcmb,yhe,annur,annunr

        h0      = dble(H_0)
        tcmb    = dble(t_cmb)
        yhe     = dble(Y_he)
        annur   = dble(a_nu_r)
        annunr  = dble(a_nu_nr)

        akmaxt=akmaxt*(h0/100.0d0)
c       if (h0.lt.25.d0.or.h0.gt.100.d0) then
c         write(*,*)
c     2      '  Warning: H0 has units of km/s/Mpc.  Your value is weird.'
c       end if
        if (tcmb.lt.2.7d0.or.tcmb.gt.2.8d0) then
          write(*,*)
     2      '  Warning: Tcmb has units of K.  Your value is weird.'
        end if

        if (yhe.lt.0.2d0.or.yhe.gt.0.3d0) then
          write(*,*)
     2      '  Warning: Y_He is the Helium fraction of baryons.',
     3      '  Your value is weird.'
        end if
        if (annunr.lt.0.or.annunr.gt.3.1) then
          write(*,*)
     2      'Warning: N_nu(massive) is strange'
          write(*,*) '  Illegal choice.  Try again:'
          go to 2
        end if
        if (annur.lt.0.or.annur.gt.3.1) then
          write(*,*)
     2      'Warning: N_nu(massless) is strange'
          write(*,*) '  Illegal choice.  Try again:'
          go to 2
        end if
        if (annunr.lt.1.and.omegan.gt.0.0) then
          write(*,*)
     2      'Warning: N_nu(massive) should be 1, 2, or 3',
     3      'For non cero omegan'
           write(*,*) '  Illegal choice.  Try again:'
          go to 2
        end if

c        write(*,*)  'Enter 0 for Peebles recombination'
c        write(*,*)  'or 1 for recfast'
c        read(*,*) rcflag
        rcflag=0
        if (rcflag.eq.1) then
           write(*,*) 'Calling recfast'
           call recfast(omegab,omegac,omegav,h0,tcmb,yhe)
           write(*,*)'recfast succesfully completed'
        endif

c       write (*,*) 'Enter 0 for no reionization'
c       write (*,*) 'Enter 1 for specified optical depth to lss (xe=1)'
c       write (*,*) 'Enter 2 for specified redshift and xe'
c       read (*,*) riflag

        if (tau.eq.0.) then
          riflag = 0
        else
          riflag = 1
        end if


c  Evaluate general cosmology constants.

c  grho gives the contribution to the expansion rate from: (g) photons,
c  (r) one flavor of relativistic neutrino (2 degrees of freedom),
c  (m) nonrelativistic matter (for Omega=1).  grho is actually
c  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
c  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
c  (Used only to set the initial conformal time.)
        grhom=3.3379d-11*h0*h0
        grhog=1.4952d-13*tcmb**4
        grhor=3.3957d-14*tcmb**4

c  adotrad gives the relation a(tau) in the radiation era:
        adotrad=2.8948d-7*tcmb*tcmb


        zri=0.0d0
        rif=0.0d0
        optdlss=0.0d0
        if (riflag.eq.1) then
c          write (*,*) 'Enter optical depth to lss'
c          read (*,*) optdlss
           optdlss = tau
           rif=1.0d0
c Calculating when reinization starts.
c           call reiopar(optdlss,zri,zristp,rif)
        end if

        if (riflag.eq.2) then
           write (*,*)'Enter redshift, ionization fraction(e.g. 50 0.2)'
           read (*,*) zri,rif
           zristp=0.07d0*zri-1.0d0
           if (zristp.lt.0.0) zristp=0.0d0
        end if

        if (ict.ne.1) then
c             write (*,*) 'Enter 0 for scalar alone,'
c             write(*,*) '1 for tensor+scalar or 2 for tensors alone'
c             read(*,*)itflag
              if (type.eq.0) then
                itflag = 0
              else
                itflag = 2
              end if

           if (itflag.ne.2) then
c             write(*,*)'number and values of scal. spec. ind. n(1,1)'
c             read(*,*)nn,(an(i),i=1,nn)
              nn=num_slope

              do i=1,num_slope
                 an(i)=slope_n(i)
                 alphans(i)=slope_al(i)
              end do
           end if
           if (itflag.eq.1) then
              write(*,*)
     2 'Tensor spectral index given by nt=ns-1 (0) or different (1)'
              read(*,*)itn
              if (itn.eq.0) then
                 do 123 in=1,nn
                    ant(in)=an(in)-1
 123             continue
              else
                 write(*,*)'values of tensor spectral indexes:'
                 read(*,*)(ant(in),in=1,nn)
              endif
           endif
           if (itflag.eq.2) then
c             write(*,*)'number and values of tens. spec. ind. (1,0)'
c             read(*,*)nn,(ant(in),in=1,nn)

              nn=num_slope
              do i=1,num_slope
                 ant(i)=slope_n(i)
              end do

              do in=1,nn
                 rat(in)=1.0d0
              enddo
           end if
           if (itflag.eq.1) then
              write(*,*)'ratio of tensor to scalar quadrupole given by'
              write(*,*)'7(1-n_S) (0) or different (1):'
              read(*,*)irt
              if (irt.eq.0) then
                 do in=1,nn
                    rat(in)=7.0d0*(1.0d0-an(in))
                 enddo
              else
                 write(*,*)'values of T/S:'
                 read(*,*)(rat(in),in=1,nn)
              endif
           end if
           if (itflag.eq.0) then
              do in=1,nn
                 rat(in)=1.0d0
              enddo
           endif

           lensflag=0
           if ((itflag.ne.2).and.(ict.eq.2)) then
              write(*,*)'Enter (0) unlensed Cls only'
              write(*,*)'Enter (1) lensed Cls, linear evolution'
              write(*,*)'Enter (2) lensed Cls, non-linear evolution'
              if (itflag.eq.1) then
                 write(*,*)'ONLY SCALAR Cls ARE LENSED'
              end if
              read(*,*)lensflag
           end if
           if ((lensflag.ne.0).and.(ict.eq.0)) then
              write(*,*)'You did not request the transfer function'
              write(*,*)'calculation needed to do the lensing'
              write(*,*)'you will have to start again'
              stop
           else
c The last requested transfer function will be used for the lensing
              ntlens=ntf
           end if

           if (itflag.ne.2) then
c             write(*,*)'Enter output filename for SCALAR cl'
              if (lensflag.ne.0) then
                 write(*,*)'If lensing was requested this will'
                 write(*,*)'be the unlensed power spectrum'
              end if

c             read(*,'(a500)')filename
              filename = outfile

              open(unit=8,file=filename,
     &             status='unknown',form='formatted')
              rewind 8

              if (lensflag.ne.0) then
                 write(*,*)'Enter output filename
     $                for LENSED SCALAR cl'
                 read(*,'(a500)')filename
                 open(unit=10,file=filename,
     &            status='unknown',form='formatted')
                 rewind 10

                 write(*,*)'Enter file with bessel functions'
                 write(*,*)' for lensing calculation (jlens.dat)'
                 read(*,'(a500)')filejlens
              end if
           end if

           if (itflag.ne.0) then
c             write(*,*)'Enter output filename for TENSOR cl'
c             read(*,'(a500)')filename
              filename = outfile
              open(unit=9,file=filename,
     &             status='unknown',form='formatted')
              rewind 9
           end if

        else
           itflag=0
           nn=0
           lensflag=0
        end if

        if (itflag.ne.2) then
           if (abs(omegak).gt.1.0d-3) then
c             write (*,*) 'CURRENTLY ONLY ISENTROPIC INITIAL CONDITIONS'
c             write(*,*) 'ARE  AVAILABLE FOR OPEN MODELS'
c             initfl=1
c             write (*,*) 'Enter initial conditions'
c             write (*,*) '1= Isentropic (adiabatic)'
c             write (*,*) '2= Isocurvature CDM'
c             write (*,*) '3= Isocurvature baryon'
c             write (*,*) '4= Isocurvature seed conditions'
c             read (*,*) initfl
              initfl = 1
           else
c             write (*,*) 'Enter initial conditions'
c             write (*,*) '1= Isentropic (adiabatic)'
c             write (*,*) '2= Isocurvature CDM'
c             write (*,*) '3= Isocurvature baryon'
c             write (*,*) '4= Isocurvature seed conditions'
c             read (*,*) initfl
              initfl = 1
           end if
        else
           initfl=1
        end if

c     TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(50,*)actual-timeprev,' END INPUT, NOW MAIN'
c       write(50,*)'PARAMETERS'
c       write(50,*)omegab,omegac,omegav,omegan
c       write(50,*)h0,tcmb,yhe,annur,annunr
c       write(50,*)riflag,zri,rif,optdlss
c       write(*,*)actual-timeprev,' END INPUT, NOW MAIN'
c       write(*,*)'PARAMETERS'
c       write(*,*)omegab,omegac,omegav,omegan
c       write(*,*)h0,tcmb,yhe,annur,annunr
c       write(*,*)riflag,zri,rif,optdlss
c     TESTING INPUT cosmological parameters.
c       write(123,*)omegab,omegac,omegav,omegan,omegak,h0,
c     &                     tcmb,yhe,annur,annunr

c call main subroutine

CC2 If you are interested in constructing a driver for
CC2 only the flat or open models code, leave the appropiate
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
           end do
        end do


       if (abs(omegak).le.1.0d-3) then
CC1 initjl reads the tables of spherical Bessel functions
CC1 so if you are making a loop over models you should
CC1 move this call outside the loop.
          if (ict.ne.1) then
c            write(*,*)'Enter input filename for jl'
c            read(*,'(a500)')filename
             filename='jl.dat'
             call initjl(filename)
          end if
          call cmbflat(clts,cltt,cles,clet,
     &    clbt,clcs,clct,clkk,cltk)
       else
             if (ict.ne.1) then
c               write(*,*)'Enter input filename for ujl'
c               read(*,'(a500)')filename
                filename='ujl.dat'
                call initujl0(filename)
             end if
             call cmbopen(clts,cltt,cles,clet,
     1     clbt,clcs,clct,clkk,cltk)
          end if

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(50,*)actual-timeprev,' END MAIN, NOW COBE'
c       write(*,*)actual-timeprev,' END MAIN, NOW COBE'


c       call COBEnormalize(clts,cltt,cles,clet,clbt,
c     &  clcs,clct,clkk,cltk)

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(50,*)actual-timeprev,' END COBE, NOW OUTPUT'


        call output(clts,cltt,cles,clet,clbt,
     &  clcs,clct,itflag,lmo)

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(50,*)actual-timeprev,' END OUTPUT, NOW LENSING'


        if (lensflag.ne.0) then
           call lensing(clts,cles,clbs,clcs,clkk,cltk)
        end if

        if (itflag.ne.2) close(8)
        if (itflag.ne.0) close(9)

c       TIMING
c       timeprev=actual
c       actual=etime(tarray)
c       write(50,*)actual-timeprev,' END LENSING, NOW END'
c       write(*,*)actual-timeinit,'sec = TOTAL TIME'

CC1 End your Loop over models here.
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
