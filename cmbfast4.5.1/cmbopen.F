C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (mzaldarriaga@cfa.harvard.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cmbopen(clts,cltt,cles,clet,clbt,
     &     clcs,clct,clkk,cltk)
c     This program integrates the linearized equations of general relativity,
c     the Boltzmann equations and the fluid equations for scalar perturbations
c     of a spatially flat Friedmann-Robertson-Walker universe in synchronous
c     gauge.  The time variable is conformal time dtau=dt/a(t) and the spatial
c     dependence is Fourier transformed.  ak is the comoving wavevector;
c     comoving distances are x=r/a(t), with a(t)=1 today.  The units of both
c     length and time are Mpc.
c
c     Presently restricted to a open (K<0) universe.
c
      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

#ifdef MPIBUILD
      include 'mpif.h'          !MPI include statements
#endif


      integer nreg(20),nr
      integer isymmfl
c
c     Dependent integration variables: a,phi,delta_c,theta_c,delta_b,theta_b,
c     2*(lmaxg+1) photon moments (total intensity plus polarization),
c     (lmaxnr+1) massless neutrino moments, (lmaxnu+1) massive neutrino
c     moments, nqmax momenta for each l.
c
c     The metric perturbation variables are computed algebraically.
      parameter (lmax0=12,lmaxnr0=25,lmaxnu0=25,nqmax0=15)
      parameter (lm2=7,lm3=4,nq1=15)
c     lmx0 MUST be grater than all of the previous lmax0.

CQ      quintessence variables

        common /qparm/ wdyn,ndyn
        double precision dynrho,wdyn_func
        external dynrho,wdyn_func,dyn_phi,dyn_nrg

      parameter (lmx0=30)
      parameter (nvar0=9+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))

c     The number will now be different
      parameter (lmaxt0=10)
      parameter (nvar0t=3*lmaxt0)
      parameter (nkmax=2*l0max)

c     Transfer functions

      double precision tautf(ntfmax),ztf(ntfmax)

c     Scalar perturbations including lensing
      double precision y(nvar0),yprime(nvar0)
      double precision d(nk0,nstep0),dpr(nk0,nstep0)
      double precision dp(nk0,nstep0),dppr(nk0,nstep0)
      double precision dkk(nk0,nstep0),dkpr(nk0,nstep0)
      double precision dl(lmax)
      double precision dpl(lmax)
      double precision dkl(lmax)
      double precision s2(nstep0,2),sp2(nstep0,2),sk2(nstep0,2)
      double precision ds2(nstep0,2),dsp2(nstep0,2),dsk2(nstep0,2)
      double precision cl(lmax,nnmax),clpr(lmax,nnmax)
      double precision cpl(lmax,nnmax),cplpr(lmax,nnmax)
      double precision ccl(lmax,nnmax),cclpr(lmax,nnmax)
      double precision ckkl(lmax,nnmax),ckklpr(lmax,nnmax)
      double precision ctkl(lmax,nnmax),ctklpr(lmax,nnmax)

c     Tensor Perturbations
      double precision yt(nvar0t),ytprime(nvar0t)
      double precision dt(nk0,nstep0),dtpr(nk0,nstep0)
      double precision dte(nk0,nstep0),dtepr(nk0,nstep0)
      double precision dtb(nk0,nstep0),dtbpr(nk0,nstep0)
      double precision dtl(lmax)
      double precision detl(lmax),dbtl(lmax)
      double precision st2(nstep0,2),ste2(nstep0,2)
      double precision stb2(nstep0,2)
      double precision dst2(nstep0,2),dste2(nstep0,2)
      double precision dstb2(nstep0,2)
      double precision ctl(lmax,nnmax),ctlpr(lmax,nnmax)
      double precision ctel(lmax,nnmax),ctelpr(lmax,nnmax)
      double precision ctbl(lmax,nnmax),ctblpr(lmax,nnmax)
      double precision ctcl(lmax,nnmax),ctclpr(lmax,nnmax)


c     Output arrays temperature: clts, cltt ; e spectra:
c     cles, clet ; b perturbation(only tensor contrubute):
c     clbt ; cross correlation: clcs, clct.

      double precision clts(l0max,nnmax),cltt(l0max,nnmax)
      double precision cles(l0max,nnmax),clet(l0max,nnmax)
      double precision clbt(l0max,nnmax)
      double precision clcs(l0max,nnmax),clct(l0max,nnmax)
      double precision clkk(l0max,nnmax),cltk(l0max,nnmax)


      double precision atau0(nstep0),dtaureg(20)
      double precision ak0(nk0)
      double precision ak1(nkmax),dak1(nkmax)
      parameter (d0hi=1.0d40,d0lo=1.0d40)
      integer l(lmax)
      double precision xl(lmax),temp,growthfactor,akgf
      common /lvalues1/ l,l0,lmo
      common /lvalues2/ akmax0
      save /lvalues1/
      save /lvalues2/
      double precision power(nk0*6,ntfmax)
      integer istore(ntfmax)
      common/transfer_function/power,istore


c     TESTING. This arrays can be used output the contributions
c     of different ks to the cls.
c      double precision ctt(lmax0),cpp(lmax0)

      common /epsilon/ epsw
c
c
c     const=7*pi**4/120.
      parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
      parameter (pio2=3.14159265d0/2.0d0,fourpi=4.0d0*3.14159265d0)
c     dverk  parameters.
      parameter (tol=1.0d-8,nw=nvar0,nwt=nvar0t)
      dimension c(24),w(nw,9)
      dimension ct(24),wt(nwt,9)
c
c
      dimension denl(lmx0),dlfdlq(nqmax0)
      common /store/ denl,dlfdlq
c
      common /timstp1/ irec,nr,nreg,itflag2
      common /timstp2/dtaureg,atau0

      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &     tcmb,yhe,annur,annunr
      common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &     nqmax,iq0,iq1,iq2
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     $     grhonr
      common /reionization/zri,taurist,zristp,tauristp,rif,optdlss

      common /lens/ lensflag
      common /curvature/ curv,r,kcurv
      double precision atrg(ntau,4),akl(nstep0),phik(nstep0)
c      common /par/atrg,taur,tau0,akl,phik,nkt

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn
      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $     itflag,lmaxt
      common /transfer/akmaxt,ztf,nlnkt,ict,ntf

      common /flip/ichiflag
      common /tnowkb/twkb1,twkb2,twkb3,twkb4

      common /maxreion/ taurmax,armax,jrmax
c     TESTING: have a flag that sets if that particular ak is outputed
c      common /akoutput/ koutflag,kfile

      double precision dtauda, dsoundda
      external oderivs,dtauda,dsoundda
      external oderivst

c     TIMING
c     Timing variables for testing purposes
c     Should be commented in released versions
c      real etime,actual,timeprev,tarray(2)
c      real time1,time2,timecl1,timecl2,timecl3,timecl4,timecl5
c      external etime

c     TO SAVE MEMORY
      common /memory1/ d,dpr,dp,dppr,dkk,dkpr
      common /memory2/dt,dtpr,dte,dtepr,dtb,dtbpr

c     K splitting
      common /cutk/ aksplit, kcutflag

#ifdef MPIBUILD
      CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )
      CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
#else
      myid=0
      numprocs=1
#endif



c     TIMING
c     Enterning CMBOPEN
c      actual=etime(tarray)
c     TESTING INPUT cosmological parameters.
c     write(123,*)omegab,omegac,omegav,omegan,omegak,h0,
c     &                     tcmb,yhe,annur,annunr
c     write(*,*)omegab,omegac,omegav,omegan,omegak,h0,
c     &                     tcmb,yhe,annur,annunr
c
      itflag2=itflag
c
c     zst is the value of z where the program stops the calculation.
c     Now set to z=0 for any model that uses this code.
      zst=0.0d0
c     if ((abs(omegab+omegac-1.0d0).lt.1.0d-3).and
c     &       .(zri.eq.0d0).and.(itflag.eq.0)) then
c     zst=10.0d0
c     else
c     zst=0.0d0
c     endif

c
      if (ict.ne.0) zst=min(zst,ztf(ntf))
c
      if (itflag.ne.2) then
         lmaxg=lmax0
      end if
#ifdef HP
            lmaxnr=lmaxnr0
#else
            lmaxnr=lm2
#endif

c     TEST, REMOVE LATER
c     nstep=nstep0

      if (itflag.ne.0) then
         lmaxt=lmaxt0
      else
         lmaxt=0
      end if
c
      do 133 j=1,lmx0
         denl(j)=1.0d0/dble(2*j+1)
 133  end do
c
c
c     Initialize neutrino mass and storage.
      if (annunr.eq.0.or.omegan.eq.0.0d0) then
         amnu=0.0d0
         nqmax=0
         lmaxnu=0
      else
c     amnu=m_nu*c**2/(k_B*T_nu0).
c     The 1.0e-18 is to prevent certain compilers (e.g. Sun's) from
c     thinking that a divide by zero will occur.
        amnu=const/(1.5d0*zeta3)*grhom/grhonr*omegan/
     &        dble((annunr+1.0d-18))
#ifdef HP
            nqmax=nqmax0
            lmaxnu=lmaxnu0
#else
            nqmax=nq1
            lmaxnu=lm3
#endif

         dq=1.0d0
         do 134 i=1,nqmax
            q=i*dq-0.5d0
            expq=exp(-q)
            dlfdlq(i)=-q/(1.0d0+expq)
 134     end do
      end if


c     Calculate number of equations
      if (itflag.ne.2) then
         iq0=11+2*lmaxg+lmaxnr
         iq1=iq0+nqmax
         iq2=iq1+nqmax
         if (ndyn.eq.0.or.ndyn.eq.3.or.ndyn.eq.4) then
          nvar=7+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1)
         else
          nvar=9+2*(lmaxg+1)+(lmaxnr+1)+nqmax*(lmaxnu+1)
         endif

      else
         nvar=0
      end if

      if (itflag.ne.0) then
         nvart=3*lmaxt
      else
         nvart=0
      end if
c     Initialize massive neutrinos.
      if ((annunr.ne.0).or.(omegan.ne.0.0d0)) then
         arel=1.d-3/amnu
         call initnu1(amnu)
      end if
c     Time today
      tau0=rombint(dtauda,1.0d-8,1.0d0,tol)
      epsw=100.0d0/tau0
c     Calculating when reinization starts.
      if (optdlss.gt.0.0d0) then
         call reiopar(optdlss,zri,zristp,rif)
      end if

c     Curvature radius
      hc=2.998d5/h0
      curv=-omegak/hc/hc
      r=1.0d0/sqrt(abs(curv))
      if (curv.lt.0.0) kcurv=-1
      if (curv.gt.0.0) kcurv=1
      if (curv.eq.0.0) kcurv=0

c     Distance to LSS
      if (((tau0/r).gt.3.1).and.(kcurv.eq.1))then
         write(*,*)'illegal parameters'
         write(*,*)'tau0/r > 3.1'
         stop
      end if
      dlss=r*sinhk(tau0/r)

      if (((tau0/r).gt.1.57).and.(kcurv.eq.1))then
         write(*,*)'WARNING: In closed models the l scale is'
         write(*,*)'compressed compared to flat models.'
         write(*,*)'This means our l sampling may not be good '
         write(*,*)'enough. Remember we only calculate some ls'
         write(*,*) 'and then we interpolate the rest. '
         write(*,*)'YOU HAVE TO TEST IF OUR SMAPLING IS GOOD ENOUGH.'
         write(*,*)'We provide lvalues_fine.inc that you'
         write(*,*)'can use. Rename lvalues_fine.inc lvalues.inc,'
         write(*,*)'recompile, rerun and compare results.'
         write(*,*)'You may not want to loose lvalues.inc, so '
         write(*,*)'save it with another name.'
         write(*,*)'jlgen and ujlgen should be rerun again.'
      end if

c     TESTING
c      write(50,*)'tau0,r'
c      write(50,*)tau0,r

c
c     Maximum and minimum k-values. dlnk is the logarithmic spacing
c     used for low k.

      akmax=akmax0/dlss
      akmin=0.15d0/dlss

c     TESTING
c      akmin=1.0d-4
c      akmax=0.22
#ifdef MAX
c     FOR PROJECT WITH MAX
c     To speed things up.
      akmax=akmax0/tau0
      akmin=0.15d0/tau0
c     TESTING
c      write(*,*)akmax,tau0
#endif
c     min k in closed universe is 3
      if (curv.gt.0.0) then
         akmin1=3.0d0/r
         akmin=max(akmin,akmin1)
      endif
      if(itflag.eq.0)then
         dlnk=0.1d0
      else
         dlnk=0.05d0
      end if
c
c
c     Timesteps during recombination (tentative, the actual
c     timestep is the minimum between this value and taurst/40,
c     where taurst is the time when recombination starts.
c     akmax is the maximum ak being calculated, 4.0/akmax
c     seems a very poor sampling, but akmax is bigger than
c     lmax/dlss, this is why this condition works.

      if (akmax.ne.0d0) then
         dtaurec=4.0d0/akmax
      else
         dtaurec=0.0d0
      end if

c     Stoping time
      taumax=rombint(dtauda,1.0d-8,1.0d0/(zst+1.0d0),tol)

c     Time at which reionization takes place
      if (zri.ne.0.0d0) then
         taurist=rombint(dtauda,1.0d-8,1.0d0/(1.0d0+zri),tol)
         tauristp=rombint(dtauda,1.0d-8,1.0d0/(1.0d0+zristp),tol)
      else
         taurist=tau0
         tauristp=tau0
      end if
c
      if (ict.eq.0) then
         taumin=0.001d0/akmax
      end if
      if (ict.eq.1) then
         taumin=0.001d0/akmaxt
      end if
      if (ict.eq.2) then
         taumin=0.001d0*min(1.0d0/akmax,1.0d0/akmaxt)
      end if
      taumin=min(taumin,0.1d0)

c
      if (amnu.ne.0.0d0) then
         taumin=min(taumin,arel/adotrad)
      end if
c     Initialize baryon temperature and ionization fractions vs. time.
c     This subroutine also fixes the timesteps where the sources are
c     saved in order to do the integration.
      call oinithermo(taumin,taumax,tau0,taurend,dtau0,nstep)
      n1=nreg(2)
      atau0min=atau0(1)
      atau0max=atau0(nstep)

c     taurend if the time recombination ends. If

      dlsso=r*sinhk((tau0-taurmax)/r)
c     Sound Horizon at LSS
      shor=rombint(dsoundda,1.0d-8,armax,tol)

c     K splitting
c     Compute cut wavevector
      aksplit=1.5d0/shor*aksplit

#ifdef MAX

c     FOR PROJECT WITH MAX
      open(unit=28,file='dlss.dat',
     &     status='unknown',form='formatted')
      rewind 28
      write(28,'(7E13.5)')dlsso,taurmax,tau0,shor,armax,r,aksplit
      close(28)

#else

c     Output dlss to shift spectra.
      if ((kcutflag.eq.1).or.(kcutflag.eq.-1)) then
c     Write dlss in the output file.
         write(8,*)'### First line has dlss for this model and kcutflag'
         write(8,*)'### Then normal l ClT ClE ClC output'
         write(8,*)dlsso,kcutflag
      end if
#endif


c     Decide if fliping of sources is necessary. Choose parameters.
      ichiflag=0
      iflip=1
      if ((irec.gt.1).and.(kcurv.eq.1)) then
         iflip=nreg(irec)
         ichiflag=1
      end if

c     TESTING
c      write(50,*)'iflip,atau0(iflip)'
c      write(50,*)iflip,atau0(iflip)

c
c     Calculating the times of the outputs of the transfer functions.
      if (ict.ne.0) then
         do 150 itf=1,ntf
            atf=1.0d0/(ztf(itf)+1.0d0)
            tautf(itf)=rombint(dtauda,1.0d-8,atf,tol)
            tautf(itf)=min(tautf(itf),atau0(nstep))

 150     enddo
      endif

c
c     Integration will only be carried out after z=10 for low
c     k, k < k10. If there is reionization the boundary will not be
c     z=10 but tauristp

      if (zst.ne.10.0d0) then
         t10=rombint(dtauda,1.0d-8,1.0d0/11.0d0,tol)
c         ak10=250.0d0/dlss
c     TESTING
         ak10=500.0d0/dlss

         t10=max(t10,atau0(nreg(nr)))
         call tau2n(t10,n10,atau0,dtaureg,nreg,nr)
      else
         t10=tau0
         ak10=akmax
         n10=nstep+1
      end if

c     TIMING
c      write(50,*)'t10,n10',t10,n10,atau0(n10)
c      write(50,*)'jrmax, nreg',jrmax
c      write(50,*)nreg
c      timeprev=actual
c      actual=etime(tarray)
c      write(50,*)actual-timeprev,' END INITIALIZE, NOW SOURCES'
c      write(*,*)actual-timeprev,' END INITIALIZE, NOW SOURCES'

c     Calculation of the CMB sources.
c
      nk=0
      if (ict.ne.1) then
c     set k values for which the sources for the anisotropy and
c     polarization will be calculated. For low values of k we
c     use a logarithmic spacing. (in the main code k is really beta)
c

c     The first ks will be set exponentially separated with dlnk0,
c     in then linearly spaced with separation dkn1 first
c     and dkn2 after.
         if (zri.ne.0.0d0) then
#ifdef HP
            dlnk0=1.0d0*dlnk
#else
            dlnk0=2.0d0*dlnk
#endif
         else
            dlnk0=5.0d0*dlnk
         end if
c
#ifdef HP
         dkn1=0.4d0/taurmax
         dkn2=0.8d0/taurmax
#else
         dkn1=0.8d0/taurmax
         dkn2=1.5d0/taurmax
#endif

c     In closed models compute the first
c     betas, from 3 to 17
         io=0
         if (curv.gt.0.0)then
            io=15
         end if

c     Find where log spacing is equal to dkn1.
c     Calculatiing the number of ks to be computed.

         nk1=int(log(dkn1/akmin/dlnk0)/dlnk0)+1
         nk1=max(nk1,1)
         akaux=akmin*exp((nk1-1.0d0)*dlnk0)

cxxx         if (akmax.gt.(1500.0d0/dlss)) then
cxxx            nk2=int((1500.0d0/dlss
cxxx     &           -akaux)/dkn1)+nk1+1
cxxx            nk=int((akmax-1500.0d0/dlss)/dkn2)+nk2+1
cxxx         else
cxxx            nk=int((akmax-akaux)/dkn1)+nk1+1
cxxx            nk2=nk+1
cxxx         end if

cxxx         akchange=1500.0d0/dlss

         akchange=5.0d0*3.14159265d0/shor

         if (akmax.gt.akchange) then
            nk2=int((akchange
     &           -akaux)/dkn1)+nk1+1
            nk=int((akmax-akchange)/dkn2)+nk2+1
         else
            nk=int((akmax-akaux)/dkn1)+nk1+1
            nk2=nk+1
         end if


c     Making sure the total number of ks does not exceed the
c     size of the allocated arrays.

             nkio=nk+io
         if ((nkio).gt.nk0) then
            write(*,*)
     2          'Sorry, the arrays were dimensioned for a maximum of'
            write(*,*)nk0, 'k modes.'
            write(*,*)'The model you requested needs',nkio
            write(*,*)'Please make the arrays bigger by making '
            write(*,*)'nk0 bigger where it appears'
            stop
         end if

         if (curv.gt.0.0) then
            do ik=1,io
               ak0(ik)=(dble(ik)+2.0d0)/r
            end do
         end if

c     TIMING
c         write(50,*)'ak sources',nk1,nk2,nk
c         write(50,*)dkn1,dkn2,akmax
c         write(*,*)'ak sources',nk1,nk2,nk

         nktemp=io+1
         do 25 ik=1,nk
            if (ik.le.nk1) then
               ak0(ik+io)=akmin*exp((ik-1.0d0)*dlnk0)
            else
               if (ik.gt.nk2) then
                  ak0(ik+io)=ak0(nk2+io)+(ik-nk2)*dkn2
               else
                  ak0(ik+io)=ak0(nk1+io)+(ik-nk1)*dkn1
               end if
            end if
 25      continue


c     Total number of ks (nk), includes added ones at the begining
c     in closed models.
         nk=nk+io


c     make ks discrete if positive curvature. Remember in this
c     part of the code ks correspond to betas in Mpc^{-1}.

         if (curv.gt.0.0) then
            aklim=3.0d0/r
            do  ik=1,nk
               if(ak0(ik).le.aklim) then
                  ak0(ik)=aklim
                  aklim=aklim+1.0d0/r
               else
                  ak0(ik)=int(ak0(ik)*r)/r
                  aklim=ak0(ik)+1.0d0/r
               endif
            end do
c     Set akmin back to previous value which will be used
c     to set ks in the time integral.
c            akmin=akmintemp
         endif

c
c     TESTING
c     To output average of the sources
c         open(unit=25,file='average_1.dat',status='unknown'
c     $        ,form='formatted')
c         rewind 25

         kkmin=1
         kkmax=nk

c     K splitting
         if ((kcutflag).eq.1) then
            do k=1,nk
               if (ak0(k).lt.(1.3d0*aksplit))kkmax=k
            end do
            kkmax=min(nk,kkmax+4)
         end if
         if ((kcutflag).eq.-1) then
            do k=1,nk
               if (ak0(k).lt.(0.3d0*aksplit))kkmin=k
            end do
            kkmin=max(1,kkmin-4)
          end if
c     TESTING
c         write(*,*)'kkmin,kkmax=',kkmin,kkmax

cxxxc BEGIN TESTING : compute one k mode
cxxx          ak=0.04
cxxx          ak2=ak*ak
cxxx          write(*,*)ak,sqrt(ak2+1.0d0/r/r)
cxxx          taustart=0.001d0/sqrt(ak2+1.0d0/r/r)
cxxx          taustart=min(taustart,0.1d0)
cxxx          taustart0=taustart
cxxx          call oinitial(y,taustart)
cxxx          tau=taustart
cxxx          ind=1
cxxxc     Begin timestep loop.
cxxx
cxxx          do j=1,nstep
cxxx             tauend=atau0(j)
cxxx             if (ak.gt.ak10.and.tauend.gt.t10.and.
cxxx     2            (ict.eq.0.or.tau.gt.tautf(ntf))) then
cxxx                dtest=0.0d0
cxxx                dptest=0.0d0
cxxx                dkktest=0.0d0
cxxx                tau=atau0(j)
cxxx             else
cxxx                tol1=tol
cxxx                call dverk(nvar,oderivs,tau,y,
cxxx     $               tauend,tol1,ind,c,nw,w)
cxxx                call oderivs(nvar,tau,y,yprime)
cxxx                call ooutput(nvar,y,yprime,j,tau0,tau,dtest,
cxxx     &               dptest,dkktest,phi)
cxxx             end if
cxxx             write(55,'(20E15.5)')tau,dtest,dptest,
cxxx     $            atau0(j),dble(j)
cxxx          end do
cxxxc     END TESTING

c     Loop over wavenumbers.
         do 30 ik=kkmin+myid,kkmax,numprocs

            ak=ak0(ik)
            ak2=ak*ak

c     TESTING. Output ks and betas where the source is calculated.
c            write(68,*)ik,ak,ak*r

c
c     Begin when wave is far outside horizon.
c     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c     of relativistic neutrinos.
            taustart=0.001d0/sqrt(ak2+1.0d0/r/r)
c     Make sure to start early in the radiation era. We are not including
c     neutrinos as sources of the tensor modes.
            taustart=min(taustart,0.1d0)
            taustart0=taustart
c     Start when massive neutrinos are strongly relativistic.
            if (amnu.ne.0.0d0) then
               arel=1.d-3/amnu
               taustart=min(taustart,arel/adotrad)
            end if

c     TESTING this variables can be used to compute
c     the sources averaged over time.
c            averaget=0.0d0
c            averagep=0.0d0

            if (itflag.ne.2) then
               call oinitial(y,taustart)
CQ      quintessence: initialize the scalar field variables (phi, phi')
CQ      in the y(nvar0-1) and y(nvar0) places. However, for this
CQ      particular k value, we only go up to y(nvar) where nvar < nvar0

           y(nvar-1) = y(nvar0-1)
           y(nvar)   = y(nvar0)

               tau=taustart
               ind=1
c     Begin timestep loop.
c
c     d contains the sources for the anisotropy and dp
c     for the polarization. t means tensor.
               if (ict.ne.0) itf=1

c     The instability reported by David Spergel ocurring at low k was
c     cured by significantly increasing the tolerance. This could lead to
c     complaints by dverk so hmin had to be lowered.

               tol1=dmin1(tol,1.0e-17*(ak/1.0e-5)**2)
               ind=2
               do j=1,24
                  c(j)=0.0d0
               end do
               c(3)=1.0d-8

              if (ict.eq.0) then
                  taulasttf=0.0d0
               else
                  taulasttf=tautf(ntf)
               end if

               do 10 j=1,nstep
                  tauend=atau0(j)
                  if (ak.gt.ak10.and.tauend.gt.t10.and.
     2                 (ict.eq.0.or.tau.gt.tautf(ntf))) then
                     d(ik,j)=0.0d0
                     dp(ik,j)=0.0d0
                     dkk(ik,j)=0.0d0
                  else
                     tol1=tol
                     call dverk(nvar,oderivs,tau,y,
     $                    tauend,tol1,ind,c,nw,w)
                     call oderivs(nvar,tau,y,yprime)
                     call ooutput(nvar,y,yprime,j,tau0,tau,d(ik,j),
     &                    dp(ik,j),dkk(ik,j),phi)

cxxxc     TESTING. Computing average sources.
cxxx                     ip1=min(j+1,nstep)
cxxx                     im1=max(1,j-1)
cxxx                     tdtau=(atau0(ip1)-atau0(im1))/dble(ip1-im1)
cxxx                     averaget=averaget+d(ik,j)*tdtau
cxxx                     averagep=averagep+dp(ik,j)*tdtau

c     TESTING. Output of sources at maximum of visibility function.
c                     if (j.eq.jrmax) then
c                        write(111,*)ak0(ik),d(ik,j),dp(ik,j)
c                     end if

c     Calculation of transfer functions.
 101                 continue
                     if (ict.ne.0.and.itf.le.ntf) then
                        if (j.lt.nstep.and.tauend.lt.tautf(itf).
     2                       and.atau0(j+1).gt.tautf(itf)) then
                           call dverk(nvar,oderivs,tau,y,tautf(itf)
     2                          ,tol1,ind,c,nw,w)
                        endif
c     output transfer functions for this k-values.
                        if (abs(tau-tautf(itf)).lt.1.0d-5) then
                           call outtransf(nvar,y,curv,ik,itf)
                           itf=itf+1
                           if (j.lt.nstep.and.itf.le.ntf.and.
     2                          atau0(j+1).gt.tautf(itf)) goto 101
                        endif
                     endif
                  end if
 10            end do
               if (lensflag.ne.0) then
                     call fderivs(nvar,tau,y,yprime)
                     call foutput(nvar,y,yprime,j,tau0,tau
     &                    ,temp,temp,temp,phi)
                     akl(ik)=ak
                     phik(ik)=phi
               endif

c     TESTING. Output averaged sources.
c               write(25,'(4E15.5)')ak0(ik),averaget,averagep,ak0(ik)*r

            end if
c
c     time loop for tensors if requested.
c
c     Tensors will only be calculated if k*tau< stsp.
c
c     stpt=50.0d0

c###  CHECK WHY THIS NUMBER IS DIFFERENT IN OPEN AND FLAT.
            stpt=100.0d0
c
            if (itflag.ne.0) then
               call oinitialt(yt,taustart0)
               tau=taustart0
               ind=1
c     Begin timestep loop.
c     dt contains the sources for the anisotropy and dtp
c     for the polarization. t means tensor.
               do 11 j=2,nstep-1
                  tauend=atau0(j)
                  if ((ak*tauend).gt.stpt) then
                     dt(ik,j)=0.0d0
                     dte(ik,j)=0.0d0
                     dtb(ik,j)=0.0d0
                  else
                     tol1=tol
                     call dverk(nvart,oderivst,tau,yt,tauend,tol1,
     &                    ind,ct,nwt,wt)
                     call oderivst(nvart,tau,yt,ytprime)
                     call ooutputt(nvart,yt,ytprime,j,tau0,tau,
     &                    dt(ik,j),dte(ik,j),dtb(ik,j))
                  end if
 11            end do
            end if
 30      end do
c     TESTING: Close output of average
c         close(25)
      end if
c
c     If transfer functions are requested, set the remaining k-values.
      if (ict.ne.0) then
         if (ict.eq.2) then
            if (akmaxt.gt.akmax) then
               nkt=nk+int((log(akmaxt)-log(akmax))*nlnkt)+1
               akdone=ak0(nk)
            else
               nkt=nk
            end if
         else
            nkt=int((log(akmaxt)-log(akmin))*nlnkt)+1
            akdone=akmin
         endif

         akgf=1.0d0
         if (akdone.gt.akgf) akgf=akdone*exp(1.5/nlnkt)

c     Loop over wavenumbers.
         do 31 ik=nk+1+myid,nkt,numprocs
            ak=akdone*exp(dble(ik-nk)/nlnkt)

c     Decrease hmin in dverk to improve accuracy at high k.
c            ind=1
            ind=2
            do j=1,24
               c(j)=0.0d0
            end do
            c(3)=1.0d-10

c     Begin when wave is far outside horizon.
c     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c     of relativistic neutrinos.
            taustart=0.001d0/sqrt(ak2+1.0d0/r/r)
c     Make sure to start early in the radiation era.
            taustart=min(taustart,0.1d0)
c     Start when massive neutrinos are strongly relativistic.
            if (amnu.ne.0.0d0) then
               arel=1.d-3/amnu
               taustart=min(taustart,arel/adotrad)
            end if
            ak2=ak*ak
            call oinitial(y,taustart)
            tau=taustart

#ifndef HISPDHIKMODES
            do 135 itf=1,ntf
               call dverk(nvar,oderivs,tau,y,tautf(itf),tol,
     2              ind,c,nw,w)
c
               call outtransf(nvar,y,curv,ik,itf)
 135        enddo
               if (lensflag.ne.0) then
                     call fderivs(nvar,tau,y,yprime)
                     call foutput(nvar,y,yprime,j,tau0,tau
     &                    ,temp,temp,temp,phi)
                     akl(ik)=ak
                     phik(ik)=phi
               endif
#else
            if(ak.le.akgf) then
               do 1352 itf=1,ntf
                  call dverk(nvar,oderivs,tau,y,tautf(itf),tol,
     2                 ind,c,nw,w)
                  call outtransf(nvar,y,0.0d0,ik,itf)
                  if (lensflag.ne.0.and.itf.eq.1) then
                    call fderivs(nvar,tau,y,yprime)
                    call foutput(nvar,y,yprime,j,tau0,tau
     &                    ,temp,temp,temp,phi1)
                  endif
                  if (lensflag.ne.0) then
                     call fderivs(nvar,tau,y,yprime)
                     call foutput(nvar,y,yprime,j,tau0,tau
     &               ,temp,temp,temp,phi)
                     akl(ik)=ak
                     phik(ik)=phi
                  endif
 1352          enddo
               closestknum=istore(1)-1
            else
               call dverk(nvar,oderivs,tau,y,tautf(1),tol,
     2              ind,c,nw,w)
               call outtransf(nvar,y,0.0d0,ik,itf)
               do 1353 itf=2,ntf
                  growthfactor=(omegac*power(6*closestknum+2,itf)
     $                 +omegab*power(6*closestknum+3,itf)
     $                 +omegan*power(6*closestknum+6,itf))
     $                 /(omegac*power(6*closestknum+2,1)
     $                 +omegab*power(6*closestknum+3,1)
     $                 +omegan*power(6*closestknum+6,1))
c     k
                  power(6*istore(itf)+1,itf)=power(6*istore(itf)+1,1)
c     all others are transfer functions
                  power(6*istore(itf)+2,itf)=
     $                 power(6*istore(itf)+2,1)*growthfactor
                  power(6*istore(itf)+3,itf)=
     $                 power(6*istore(itf)+3,1)*growthfactor
                  power(6*istore(itf)+4,itf)=
     $                 power(6*istore(itf)+4,1)*growthfactor
                  power(6*istore(itf)+5,itf)=
     $                 power(6*istore(itf)+5,1)*growthfactor
                  power(6*istore(itf)+6,itf)=
     $                 power(6*istore(itf)+6,1)*growthfactor
                  istore(itf) = istore(itf) +1
 1353          enddo
            endif
#endif

 31      end do
      endif

      if (lensflag.ne.0) then
c run k=akgf mode again to get finely sampled growth factor
            ak=akgf

c     Decrease hmin in dverk to improve accuracy at high k.
c            ind=1
            ind=2
            do j=1,24
               c(j)=0.0d0
            end do
            c(3)=1.0d-10

c     Begin when wave is far outside horizon.
c     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c     of relativistic neutrinos.
            taustart=0.001d0/sqrt(ak2+1.0d0/r/r)
c     Make sure to start early in the radiation era.
            taustart=min(taustart,0.1d0)
c     Start when massive neutrinos are strongly relativistic.
            if (amnu.ne.0.0d0) then
               arel=1.d-3/amnu
               taustart=min(taustart,arel/adotrad)
            end if
            ak2=ak*ak
            call oinitial(y,taustart)
            tau=taustart
            a0=1.0d-8
c change ar below for nonstandard reconbination
            ar=1.0d0/1089.0d0
            taur=rombint(dtauda,a0,ar,tol)

            do i=1,ntau
              taunew=dble(i)*(tau0-taur)/dble(ntau)+taur
              call dverk(nvar,oderivs,tau,y,taunew,tol,
     2              ind,c,nw,w)
              call oderivs(nvar,tau,y,yprime)
              call ooutput(nvar,y,yprime,j,tau0,tau
     &             ,d(ik,j),dp(ik,j),dkk(ik,j),phi)
              tau=taunew
c              rchi=tau0-tau
c              rr=tau0-taur
              atrg(i,1)=tau
              atrg(i,2)=y(1)
              atrg(i,3)=phi

            enddo
        endif





#ifdef MPIBUILD
c
c       communicate values of sources to outher processors
c
      if(myid.eq.0) write(*,*) 'starting communication'
      CALL MPI_BARRIER_(MPI_COMM_WORLD,IERR)
      CALL COMMUNICATE(d,dp,dkk,dt,dte,dtb,nstep,nk,itflag)
      CALL MPI_BARRIER_(MPI_COMM_WORLD,IERR)
#endif

c
      if(ict.ne.0) then
        call output_power(ntf,amnu)

#ifdef LYAOUTPUTS

         closestk=1e10
         closestknum=1

c     here we find k closest to k=1

         do itf = 1,istore(ntf)
            temp=abs(log(power(6*(itf-1)+1,ntf)/1.0d10))
            if (temp.ne.0) then
               if (temp.le.closestk) then
                  closestknum=itf
                  closestk=temp
               endif
            endif
         end do

         do itf=ntf,1,-1
            atf=1.0d0/(ztf(itf)+1.0d0)
            dtaudatf=dtauda(atf)
            htf=1.0d0/dtaudatf/atf**2*2.998d5
            if (ndyn.eq.0) then
               omegaq = omegav*atf**4*grhom*dtaudatf**2/3
            else
               omegaq=omegav*dynrho(atf)*atf**4*grhom*dtaudatf**2/3
            endif

c            k(ncount) = power(6*(closestknum-1)+1,itf)
c            tfc(ncount) = power(6*(closestknum-1)+2,itf)
c            tfb(ncount) = power(6*(closestknum-1)+3,itf)
c            tfg(ncount) = power(6*(closestknum-1)+4,itf)
c            tfn(ncount) = power(6*(closestknum-1)+5,itf)
c            tfnm(ncount) = power(6*(closestknum-1)+6,itf)


            write(98,'(6E14.6)')ztf(itf),
     $           (omegac*power(6*(closestknum-1)+2,itf)
     $           +omegab*power(6*(closestknum-1)+3,itf)
     $           +omegan*power(6*(closestknum-1)+6,itf))
     $           /(omegac*power(6*(closestknum-1)+2,ntf)
     $           +omegab*power(6*(closestknum-1)+3,ntf)
     $           +omegan*power(6*(closestknum-1)+6,ntf)),
     $           htf,tautf(itf),omegaq,wdyn_func(atf)
         enddo
         close(98)
#endif

      endif



c     TIMING
c      timeprev=actual
c      actual=etime(tarray)
c      write(50,*)actual-timeprev,' END SOURCES, NOW INTERPOL'
c      write(*,*)actual-timeprev,' END SOURCES, NOW INTERPOL'
c
c     if CMB calculations are requested, calculate the Cl by
c     integrating the sources over time and over k.
c
c
      if (ict.ne.1) then
         if (itflag.ne.2) then
c     get the interpolation matrix for the sources to interpolate them
c     for other k-values, scalar case.
            do 100 i=1,nstep
               call spline(ak0,d(1,i),nk,d0lo,d0hi,dpr(1,i))
               call spline(ak0,dp(1,i),nk,d0lo,d0hi,dppr(1,i))
               call spline(ak0,dkk(1,i),nk,d0lo,d0hi,dkpr(1,i))
 100        continue
            do 105 in=1,nn
               do 106 j=1,l0
                  cl(j,in)=0.0d0
                  cpl(j,in)=0.0d0
                  ccl(j,in)=0.0d0
                  ckkl(j,in)=0.0d0
                  ctkl(j,in)=0.0d0
 106           continue
 105        continue
         end if
c
c     get the interpolation matrix for the tensor sources.
         if (itflag.ne.0) then
            do 110 i=1,nstep
               call spline(ak0,dt(1,i),nk,d0lo,d0hi,dtpr(1,i))
               call spline(ak0,dte(1,i),nk,d0lo,d0hi,dtepr(1,i))
               call spline(ak0,dtb(1,i),nk,d0lo,d0hi,dtbpr(1,i))
 110        continue
            do 115 in=1,nn
               do 116 j=1,l0
                  ctl(j,in)=0.0d0
                  ctel(j,in)=0.0d0
                  ctbl(j,in)=0.0d0
                  ctcl(j,in)=0.0d0
 116           continue
 115        continue
         end if
c
c     Fixing the # of k for the integration.
c     We will use logarithmic spacing first, then dk0 and then dk.
c     We will switch from logarithmic to linear when k=10*dk0 which
c     corresponds to l=15. This is to catch the plateau at low l.
c     Then I will increase the step at k=500*dk0 which corrsponds
c     to l=750.

c     dk=2.5d0/dlss
c     dk0=1.5d0/dlss
c     no=500

#ifdef HP
         dk=1.1d0/tau0
         dk0=1.1d0/tau0
         no=int(700.0d0*tau0/dlss)
         dlnk1=0.07d0
#else
         dk=2.5d0/tau0
         dk0=1.5d0/tau0
         no=int(500.0d0*tau0/dlss)
         dlnk1=0.1d0
#endif

cxxx         write(*,*)'Enter accfact'
cxxx         read(*,*)accfact
cxxx         dk=dk/accfact
cxxx         dk0=dk0/accfact
cxxx         no=no*accfact
cxxx         dlnk1=dlnk1/accfact

c     TESTING
c         write(50,*)'dks',dk0*r,dk*r

c     In closed models making the steps in k intenger.
         if (kcurv.eq.1) then
            kd=int(dk0*r)
            kd=max(kd,1)

            dk0=dble(kd)/r

            kd2=int(dk*r)
            kd=max(kd2,1)

c     TESTING
c            kd=1

            dk=dble(kd)/r
         end if


c     Calculating the number of ks for the time integration.

         no1=int(log(10.0d0*dk0/akmin)/dlnk1)+1
         if (akmax.gt.(no*dk0)) then
            nko=int((akmax-no*dk0)/dk)+no
         else
            no=int((akmax-10.0d0*dk0)/dk0)+no1
            nko=no
         end if

c     Making sure not to exceed the array sizes.

         if (nko.gt.nkmax) then
            write(*,*)
     2          'Sorry, the arrays were dimensioned for a maximum of'
            write(*,*) nkmax, 'k modes.'
            write(*,*)'The model you requested needs',nko
            write(*,*)'Please make the arrays bigger by making '
            write(*,*)'nkmax bigger where it appears'
            stop
         end if

c     TIMING
c         write(50,*)'ak integral',no1,no,nko
c         write(50,*)dk0,dk,akmax
c         write(*,*)'ak integral',no1,no,nko

c     ks where the integral is calculated.
         do 198 k=1,nko
            if (k.le.no) then
               if (k.le.no1) then
                  ak1(k)=10.0d0*dk0*exp(-(no1-k)*dlnk1)
               else
                  ak1(k)=ak1(no1)+(k-no1)*dk0
               end if
            else
               ak1(k)=ak1(no)+(k-no)*dk
            end if
 198     continue
         dak1(1)=0.5d0*ak1(1)

c     Discrete if closed universe.
         if (curv.gt.0.0) then
            aklim=3.0d0/r
            do k=1,nko
               if(ak1(k).le.aklim) then
                  ak1(k)=aklim
                  aklim=aklim+1.0d0/r
               else
                  ak1(k)=int(ak1(k)*r+0.1d0)/r
                  aklim=ak1(k)+1.0d0/r
               endif
            end do
            dak1(1)=1.0d0/r
         endif

c     Calculating dk for the k integral.
         do k=2,nko-1
            dak1(k)=0.5d0*(ak1(k+1)-ak1(k-1))
         enddo
         dak1(nko)=0.5d0*(ak1(nko)-ak1(nko-1))

c     TESTING Output ks where the integral is calculated.
c         do k=2,nko-1
c            akt=ak1(k)
c            write(445,'(1I6,5E15.5)')k,akt,akt*r
c     $           ,dak1(k)*r,dak1(k)
c         enddo


c     klo and khi bracket the position of this particular
c     k in the k table of the sources.
         klo=1
         khi=2

c     TIMING
c         timeprev=actual
c         actual=etime(tarray)
c         write(50,*)actual-timeprev,' END INTERPOL, NOW Cls'
c         write(*,*)actual-timeprev,' END INTERPOL, NOW Cls'
c         timecl1=0.0
c         timecl2=0.0
c         timecl3=0.0
c         timecl4=0.0
c         timecl5=0.0

c     TESTING
c     To output average of the sources
c         open(unit=24,file='average_2.dat',status='unknown'
c     $        ,form='formatted')
c         rewind 24


         kkmin=1
         kkmax=nko

c     K splitting
         if ((kcutflag).eq.1) then
            do k=1,nko
               if (ak1(k).lt.(1.3d0*aksplit))kkmax=k+1
            end do
c            ak1(kkmax)=aksplit
c            dak1(kkmax)=0.5d0*(ak1(kkmax)-ak1(kkmax-1))
         end if
         if ((kcutflag).eq.-1) then
            do k=1,nko
               if (ak1(k).le.(0.3d0*aksplit))kkmin=k
            end do
c            ak1(kkmin)=aksplit
c            dak1(kkmin)=0.5d0*(ak1(kkmin+1)-ak1(kkmin))
         end if
c         write(*,*)'kkmin,kkmax=',kkmin,kkmax

c     Begin k-loop
         do 200 k=kkmin+myid,kkmax,numprocs
            akt=ak1(k)

c     TESTING : DECIDE if I will output the ujl as a function of time for
c     this beta
c            koutflag=0
c            if(k.eq.1) then
c               koutflag=1
c               kfile=460
c            end if
c            if(k.eq.2) then
c               koutflag=1
c               kfile=461
c            end if
c            if(k.eq.3) then
c               koutflag=1
c               kfile=462
c            end if
c            if(k.eq.4) then
c               koutflag=1
c               kfile=463
c            end if
c            if(k.eq.5) then
c               koutflag=1
c               kfile=464
c            end if
c            if(k.eq.6) then
c               koutflag=1
c               kfile=465
c            end if
c            if(k.eq.954) then
c               koutflag=1
c               kfile=466
c            end if
c            if(k.eq.955) then
c               koutflag=1
c               kfile=467
c            end if
c            if(k.eq.956) then
c               koutflag=1
c               kfile=468
c            end if
c            if(k.eq.957) then
c               koutflag=1
c               kfile=469
c            end if
c            if(k.eq.958) then
c               koutflag=1
c               kfile=470
c            end if


            beta=akt*r
            betam1=1.0d0/beta



c     finding position of k in table ak0 to do the interpolation.
 211        continue
            if ((akt.gt.ak0(klo+1)).and.(klo.lt.(nk-1))) then
               klo=klo+1
               khi=klo+1
               goto 211
            end if
            ho=ak0(khi)-ak0(klo)
            a0=(ak0(khi)-akt)/ho
            b0=(akt-ak0(klo))/ho
            if(itflag.ne.2)  then
               do 250 j=1,l0
                  dl(j)=0.0d0
                  dpl(j)=0.0d0
                  dkl(j)=0.0d0
 250           end do

c     TIMING
c               time1=etime(tarray)

c
c     Interpolating the source as a function of time for the present
c     wavelength.
               if (akt.lt.ak10) then
                  nstps=nstep-1
               else
                  nstps=n10
               end if

c     TESTING. Compute averaged sources for each k in the time
c     integral to compare with those calculated earlier.
c               averaget=0.0d0
c               averagep=0.0d0

               do 304 i=1,nstps
                  s2(i,1)=a0*d(klo,i)+b0*d(khi,i)
     $                 +((a0**3-a0)*dpr(klo,i)
     &                 +(b0**3-b0)*dpr(khi,i))*ho*ho/6.d0
                  sp2(i,1)=a0*dp(klo,i)+b0*dp(khi,i)+((a0**3-a0)
     &                 *dppr(klo,i)
     $                 +(b0**3-b0)*dppr(khi,i))*ho*ho/6.d0
                  sk2(i,1)=a0*dkk(klo,i)+b0*dkk(khi,i)+((a0**3-a0)
     &                 *dkpr(klo,i)
     $                 +(b0**3-b0)*dkpr(khi,i))*ho*ho/6.d0

c     TESTING. Compute averaged sources for each k in the time
c     integral to compare with those calculated earlier.
c                  ip1=min(i+1,nstps)
c                  im1=max(1,i-1)
c                  tdtau=(atau0(ip1)-atau0(im1))/dble(ip1-im1)
c                  averaget=averaget+s2(i,1)*tdtau
c                  averagep=averagep+sp2(i,1)*tdtau
c                  if (koutflag.eq.1) then
c                     write(kfile+100,'(3E15.5,1I6)')atau0(i),s2(i,1),
c     $                    (tau0-atau0(i))/r,i
c                  end if

 304           continue
               do i=nstps+1,nstep
                  s2(i,1)=0.0d0
                  sp2(i,1)=0.0d0
                  sk2(i,1)=0.0d0
               end do
c     TESTING
c            if (koutflag.eq.1) then
c               do i=1,nstep
c                  write(kfile+300,'(3E15.5)')atau0(i),s2(i,1),
c     $                 (atau0(nstep)-atau0(i))/r
c               end do
c            end if


c     If need to symmetrize sources because chi exceeds pi/2
               if (ichiflag.eq.1) then
                  call symsource(s2,iflip,nstep)
                  call symsource(sp2,iflip,nstep)
                  call symsource(sk2,iflip,nstep)
               end if


c     TESTING. Output averaged sources.
c               averaget2=0.0d0
c               averagep2=0.0d0
c               averaget3=0.0d0
c               averagep3=0.0d0
c               do i=iflip,nstps
c                  ip1=min(i+1,nstps)
c                  im1=max(iflip,i-1)
c                  tdtau=(atau0(ip1)-atau0(im1))/dble(ip1-im1)
c                  averaget2=averaget2+s2(i,1)*tdtau
c                  averagep2=averagep2+sp2(i,1)*tdtau
c                  averaget3=averaget3+s2(i,2)*tdtau
c                  averagep3=averagep3+sp2(i,2)*tdtau
c               end do
c               write(24,'(8E15.5)')akt
c     $              ,averaget,averagep
c     $              ,averaget2,averagep2
c     $              ,averaget3,averagep3,beta

c     TESTING
c            if (koutflag.eq.1) then
c               do i=1,nstep
c                  write(kfile+200,'(4E15.5)')atau0(i),
c     $                 s2(i,1),s2(i,2),
c     $                 (atau0(nstep)-atau0(i))/r
c               end do
c            end if
c

c     Getting interpolation table for the source at this time.
               call spline(atau0,s2(1,1),nstep,d0lo,d0hi,ds2(1,1))
               call spline(atau0,sp2(1,1),nstep,d0lo,d0hi,dsp2(1,1))
               call spline(atau0,sk2(1,1),nstep,d0lo,d0hi,dsk2(1,1))

c     Also interpolate antisymmetric sources if necessary
               if (ichiflag.eq.1) then
                  call spline(atau0,s2(1,2),nstep,
     $                 d0lo,d0hi,ds2(1,2))
                  call spline(atau0,sp2(1,2),nstep,
     $                 d0lo,d0hi,dsp2(1,2))
                  call spline(atau0,sk2(1,2),nstep,
     $                 d0lo,d0hi,dsk2(1,2))
               end if

            end if

c     TIMING
c            time2=etime(tarray)
c            timecl1=time2-time1+timecl1

c     If tensors wanted
            if (itflag.ne.0) then
               do 260 j=1,l0
                  dtl(j)=0.0d0
                  detl(j)=0.0d0
                  dbtl(j)=0.0d0
 260           end do


c
c     Interpolating the tensor source as a function of time for the present
c     wavelength.
               st2(1,1)=0.0d0
               ste2(1,1)=0.0d0
               stb2(1,1)=0.0d0
               st2(nstep-1,1)=0.0d0
               ste2(nstep-1,1)=0.0d0
               stb2(nstep-1,1)=0.0d0
               st2(nstep,1)=0.0d0
               ste2(nstep,1)=0.0d0
               stb2(nstep,1)=0.0d0
               nspt=2
               do 306 i=2,nstep
                  xf=akt*(tau0-atau0(i))

                  if (((akt*atau0(i)).lt.stpt).and.
     $                 (xf.gt.1.0d-8)) then
                     nspt=i
                     st2(i,1)=a0*dt(klo,i)+b0*dt(khi,i)+((a0**3-a0)
     &                    *dtpr(klo,i)+(b0**3-b0)
     $                    *dtpr(khi,i))*ho*ho/6.d0
                     ste2(i,1)=a0*dte(klo,i)+b0*dte(khi,i)
     $                    +((a0**3-a0)
     &                    *dtepr(klo,i)+(b0**3-b0)
     $                    *dtepr(khi,i))*ho*ho/6.d0
                     stb2(i,1)=a0*dtb(klo,i)+b0*dtb(khi,i)
     $                    +((a0**3-a0)
     &                    *dtbpr(klo,i)+(b0**3-b0)
     $                    *dtbpr(khi,i))*ho*ho/6.d0
                  else
                     st2(i,1)=0.0d0
                     ste2(i,1)=0.0d0
                     stb2(i,1)=0.0d0
                  end if
 306           end do
               nspt=max(nspt,n1)

c     If need to symmetrize sources because chi exceeds pi/2
               if (ichiflag.eq.1) then
                  call symsource(st2,iflip,nstep)
                  call symsource(ste2,iflip,nstep)
                  call symsource(stb2,iflip,nstep)
               end if

c     Getting interpolation table for the source at this time.
               call spline(atau0,st2(1,1),nstep,d0lo,d0hi,dst2(1,1))
               call spline(atau0,ste2(1,1),nstep,d0lo
     $              ,d0hi,dste2(1,1))
               call spline(atau0,stb2(1,1),nstep,d0lo
     $              ,d0hi,dstb2(1,1))

c     Also interpolate antisymmetric sources if necessary
               if (ichiflag.eq.1) then
                  call spline(atau0,st2(1,2),nstep,
     $                 d0lo,d0hi,dst2(1,2))
                  call spline(atau0,ste2(1,2),nstep,
     $                 d0lo,d0hi,dste2(1,2))
                  call spline(atau0,sk2(1,2),nstep,
     $                 d0lo,d0hi,dstb2(1,2))
               end if


            end if
c
c     Begin l and  time-loop to integrate scalar and tensor perturbations.
c
c     Determining ranges of integration


            do 310 j=1,l0
c     Setting initial conditions for calculating ujl.

               ll=l(j)

c     For closed models beta has to be larger than l, is not the
c     ujl is zero. If beta is smaller or equal I will skip the
c     rest of the ls. To avoid problems with int, which sometimes
c     does funny things when beta is intenger.0d0 add 0.1 to beta
c     before computing int.

               if ((ll.ge.int(beta+0.1d0)).and.(kcurv.eq.1))goto 400

c     Decide if we are going to use symmetric or anti-symmetric
c     source.
c     If beta-l-1 is even the function is symmetric else is
c     it is antisymmetric
               isymmfl=1
               iphase=1
               if (ichiflag.eq.1) then
                  iphase=0
                  ioddeven=mod(int(beta+0.1d0)-ll,2)
                  if (ioddeven.eq.0) isymmfl=2
               end if


c     TESTING
c               if (ll.eq.950) then
c                  write(231,*)int(beta+0.1),ioddeven,isymmfl
c               end if

c     For this beta and l find where to start the ujl integration,
c     find the values of y1=ujl*r and its derivative y2.

               call initchi(beta,j,ll,chi0,y1,y2,kcurv)

c       Changed to fix bug
                y1t=y1
                y2t=y2



c     finding the time ranges  for this beta,l from the table

               tmax=(tau0-r*chi0)
c     TESTING
c               write(50,*)'tmax,tau0,r,chi0'
c               write(50,*)tmax,tau0,r,chi0

               if (tmax.lt.atau0min)goto 400
               tmax=min(tmax,atau0max)

               call tau2n(tmax,nend,atau0,dtaureg,nreg,nr)
               nend=min(nend,nstep-5)


c     Finding the place where we will approximate by wkb formula.

               aux1=sinhk(chi0)
c     TESTING
c               aux1=2.0*aux1
               if ((ll.gt.10).and.(ll.lt.400)) then
                  aux1=8.0d0*aux1
               else
                  aux1=3.0d0*aux1
               end if
               if (kcurv.eq.-1) then
                  chimax=log(aux1+sqrt(aux1**2+1.0d0))
                  chimax=max(chi0+100.0d0*betam1,chimax)
               else
                  aux1=min(aux1,1.0d0)
                  chimax=asin(aux1)
                  chimax=max(chi0+100.0d0*betam1,chimax)
c     To be safer, although have not found a model where it mattered,
c     if the flipping ocurred delay WKB approximation. In any case
c     if you  want a crazy models you should pay the price in performance.
                  if (ichiflag.eq.1) chimax=max(chi0
     $                 +300.0d0*betam1,chimax)

c     TESTING
c                  if (ichiflag.eq.1) chimax=max(chi0
c     $                 +900.0d0*betam1,chimax)

               end if

               tmin1=(tau0-r*chimax)
               if ((tmin1.gt.twkb1).and.(tmin1.lt.twkb2))tmin1=twkb1
               if ((tmin1.gt.twkb3).and.(tmin1.lt.twkb4))tmin1=twkb3

               tmin1=max(tmin1,atau0min)
               tmin1=min(tmin1,atau0max)

c     TESTING
c               write(50,*)'tmin1,tau0,r,chimax'
c               write(50,*)tmin1,tau0,r,chimax


               call tau2n(tmin1,nst1,atau0,dtaureg,nreg,nr)
               nst1=min(nst1,nend-1)
               nst1=max(nst1,1)
c     TESTING
c               if (ll.eq.30) then
c                  write(232,'(4E15.5,1I6)')akt*r
c     $                 ,chi0+300.0d0*betam1,chimax,asin(aux1),nst1
c               end if

c     TESTING
c               if (koutflag.eq.1) then
c                  write(419,*)k,ll,nst1,nend,0
c               end if
c     If for a particular l and beta we have to start before recombination
c     for the Bessel function to be different from zero, the situation
c     will be even worse for higher l, so skip higher ls.

               if (nend.lt.2) goto 400


c     Scalar modes

               if (itflag.ne.2) then

                  ntmp=max(nstps,2*iflip)
                  if (nst1.ge.ntmp) then
                     nst1=nstps
                     ndone=nstps

                     y1=0
                     y2=1.0d0
c     Skip the range where ujl is being solved for and
c     just use the wkb approximation
                     goto 280
                  end if
c     TIMING
c                  time1=etime(tarray)

c     Calculate initial conditions.
                  call initujl(ll,akt,chi0,atau0(nend),y1,
     &                 y2,r,tau0)
c
                  ndone=nend
c     TIMING
c                  time2=etime(tarray)
c                  timecl2=time2-time1+timecl2

c     Integration when ujl is being calculated by integrating the
c     differential equation.

c     TIMING
c                  time1=etime(tarray)

                  do ir=nr,irec,-1
                     nb=nreg(ir+1)
                     na=nreg(ir)
c     TESTING
c                     if (koutflag.eq.1) then
c                        write(419,*)k,ll,na,nb,1
c                     end if

                     if ((na.lt.nend).and.(nb.gt.nst1)) then
                        nst2=max(nst1,na)
c     TESTING
c                        if (koutflag.eq.1) then
c                           write(419,*)k,ll,ndone,nst2,2
c                        end if

                        call intopen1(atau0,nst2,ndone,dtaureg(ir)
     1                       ,akt,ll,tau0,r,s2(1,isymmfl)
     $                       ,sp2(1,isymmfl),sk2(1,isymmfl),
     2                       ds2(1,isymmfl),dsp2(1,isymmfl),
     $                       dsk2(1,isymmfl),y1,y2,out1,out2,out3)
                        dl(j)=dl(j)+out1
                        dpl(j)=dpl(j)+out2
                        dkl(j)=dkl(j)+out3
                        ndone=nst2
                     end if
                  end do

c     TIMING
c                  time2=etime(tarray)
c                  timecl3=time2-time1+timecl3


 280              continue

c     Integration when ujl is being approximated by its asimptotic form.
c     differential equation.

c     TIMING
c                  time1=etime(tarray)

                  do ir=nr,irec,-1
                     nb=nreg(ir+1)
                     na=nreg(ir)
                     if (nst1.gt.na) then
c     TESTING
c                     if (koutflag.eq.1) then
c                        write(419,*)k,ll,na,nb,3
c                     end if


                        call intopen2(atau0,na,ndone,dtaureg(ir),
     1                       akt,ll,tau0,r,s2(1,isymmfl),sp2(1,
     $                       isymmfl),sk2(1,isymmfl)
     $                       ,ds2(1,isymmfl),dsp2(1,isymmfl)
     $                       ,dsk2(1,isymmfl),
     2                       y1,y2,out1,out2,out3)
                        dl(j)=dl(j)+out1
                        dpl(j)=dpl(j)+out2
                        dkl(j)=dkl(j)+out3
                        ndone=na
                     end if
                  end do
               end if



c     TIMING
c               time2=etime(tarray)
c               timecl4=time2-time1+timecl4

c     Tensor modes

               if (itflag.ne.0) then
                  ntmp=max(nspt,2*iflip)
                  y1=y1t
                  y2=y2t
                  if (nst1.ge.ntmp) then
                     nst1=nspt
                     ndone=nspt

                     y1=0
                     y2=1.0d0
c     skip the part where the ujl is being solved for and
c     use only wkb approximation.
                     goto 290
                  end if

c     Calculate initial conditions.
                  call initujl(ll,akt,chi0,atau0(nend),y1,
     &                 y2,r,tau0)
c
                  ndone=nend

c     Integration when ujl is being calculated by integrating the
c     differential equation.

                  do ir=nr,irec,-1
                     nb=nreg(ir+1)
                     na=nreg(ir)

                     if ((na.lt.nend).and.(nb.gt.nst1)) then
                        nst2=max(nst1,na)
                        call intopen1t(atau0,nst2,ndone
     &                       ,dtaureg(ir),akt,ll,tau0
     $                       ,r,st2(1,isymmfl)
     $                       ,ste2(1,isymmfl),stb2(1,isymmfl)
     $                       ,dst2(1,isymmfl)
     $                       ,dste2(1,isymmfl),dstb2(1,isymmfl),
     &                       y1,y2,out1,out2,out3)
                        dtl(j)=dtl(j)+out1
                        detl(j)=detl(j)+out2
                        dbtl(j)=dbtl(j)+out3
                        ndone=nst2
                     end if
                  end do


 290              continue

c     Integration when ujl is being approximated by its asimptotic form.
c     differential equation.

                  do ir=nr,irec,-1
                     nb=nreg(ir+1)
                     na=nreg(ir)

                     if (nst1.gt.na) then
                        call intopen2t(atau0,na,ndone
     &                       ,dtaureg(ir),akt,ll,tau0,r
     $                       ,st2(1,isymmfl)
     &                       ,ste2(1,isymmfl),stb2(1,isymmfl)
     $                       ,dst2(1,isymmfl)
     $                       ,dste2(1,isymmfl),dstb2(1,isymmfl),
     &                       y1,y2,out1,out2,out3)
                        dtl(j)=dtl(j)+out1
                        detl(j)=detl(j)+out2
                        dbtl(j)=dbtl(j)+out3
                        ndone=na
                     end if
                  end do
               end if


c     This is the end of the loop over l.
 310        end do

 400        continue

c     TESTING. Output the result of line of sight integral
c     for several ls as a function of k.
c            write(350,'(5E15.5)')akt,dl(1),dl(2),dl(3),dl(4)
c            write(351,'(5E15.5)')akt,dl(5),dl(6),dl(7),dl(8)
c            write(352,'(5E15.5)')akt,dl(9),dl(10),dl(11),dl(12)
c            write(353,'(5E15.5)')akt,dl(13),dl(14),dl(15),dl(16)
c            write(354,'(5E15.5)')akt,dl(17),dl(18),dl(19),dl(20)
c            write(355,'(5E15.5)')akt,dl(21),dl(22),dl(23),dl(24)
c            write(356,'(5E15.5)')akt,dl(25),dl(26),dl(27),dl(28)
c            write(357,'(5E15.5)')akt,dl(29),dl(30),dl(31),dl(32)
c            write(358,'(5E15.5)')akt,dl(33),dl(34),dl(35),dl(36)
c            write(359,'(5E15.5)')akt,dl(37),dl(38),dl(39),dl(40)

c     Adding to calculate the integral over k.
c     Scalar case

c
c     TIMING
c            time1=etime(tarray)

            if (itflag.ne.2) then
               do in=1,nn
                  call powersopen(akt,in,apowers)
                  apowers=apowers/akt
                  do j=1,l0
                     ckj=apowers*dl(j)*dl(j)*dak1(k)
                     cpkj=apowers*dpl(j)*dpl(j)*dak1(k)
                     cckj=apowers*dl(j)*dpl(j)*dak1(k)
                     ckkkj=apowers*dkl(j)*dkl(j)*dak1(k)
                     ctkkj=apowers*dl(j)*dkl(j)*dak1(k)

                     cl(j,in)=cl(j,in)+ckj
                     cpl(j,in)=cpl(j,in)+cpkj
                     ccl(j,in)=ccl(j,in)+cckj
                     ckkl(j,in)=ckkl(j,in)+ckkkj
                     ctkl(j,in)=ctkl(j,in)+ctkkj

c     TESTING. Keep integrand of cls for each k.
c                     ctt(j)=ckj
c                     cpp(j)=cpkj
                  end do
               end do

c     TESTING. Output Integrand of cl as a function of k for
c     several ls.
c               write(450,'(5E15.5)')akt,ctt(1),ctt(2),ctt(3),ctt(4)
c               write(451,'(5E15.5)')akt,ctt(5),ctt(6),ctt(7),ctt(8)
c               write(452,'(5E15.5)')akt,ctt(9),ctt(10)
c     $              ,ctt(11),ctt(12)
c               write(453,'(5E15.5)')akt,ctt(13)
c     $              ,ctt(14),ctt(15),ctt(16)
c               write(454,'(5E15.5)')akt,ctt(17)
c     $              ,ctt(18),ctt(19),ctt(20)
c               write(455,'(5E15.5)')akt,ctt(21),ctt(22),ctt(23),ctt(24)
c               write(456,'(5E15.5)')akt,ctt(35),ctt(26),ctt(27),ctt(28)
c               write(457,'(5E15.5)')akt,ctt(29),ctt(30),ctt(31),ct(32)
c               write(458,'(5E15.5)')akt,ctt(33),ct(34)
c     $              ,ctt(35),ctt(36)
c               write(459,'(5E15.5)')akt,ctt(37),ctt(38)
c     $              ,ctt(39),ctt(40)
c
            end if

c     TIMING
c            time2=etime(tarray)
c            timecl5=time2-time1+timecl5

c     Tensor case
c

            if (itflag.ne.0) then
               do in=1,nn
                  call powertopen(akt,in,apowert)
                  apowert=apowert/akt
                  do j=1,l0
                     ctkj=apowert*dtl(j)*dtl(j)*dak1(k)
                     ctl(j,in)=ctl(j,in)+ctkj

                     ctekj=apowert*detl(j)*detl(j)*dak1(k)
                     ctel(j,in)=ctel(j,in)+ctekj

                     ctbkj=apowert*dbtl(j)*dbtl(j)*dak1(k)
                     ctbl(j,in)=ctbl(j,in)+ctbkj

                     ctckj=apowert*dtl(j)*detl(j)*dak1(k)
                     ctcl(j,in)=ctcl(j,in)+ctckj
                  end do
               end do
            end if

c     This is the end of the loop over beta!
 200     end do
#ifdef MPIBUILD
      CALL TRANSFER_CL(cl,cpl,ccl,ckkl,ctkl,ctl,ctel,ctbl,ctcl,
     >                  itflag,nn,l0)
        IF(MYID.ne.0) return
#endif

c     TESTING: Close file for output average
c         close(24)

c     TIMING
c         timeprev=actual
c         actual=etime(tarray)
c         write(50,*)actual-timeprev,' END Cls, NOW FINAL'
c         write(50,*)timecl1,timecl2,timecl3
c         write(50,*)timecl4,timecl5
c         write(*,*)actual-timeprev,' END Cls, NOW FINAL'
c         write(*,*)timecl1,timecl2,timecl3
c         write(*,*)timecl4,timecl5

c
c
c     Final calculations for CMB output.
c     Normalize the spectra and interpolate to have every l.

         do j=1,l0
            xl(j)=dble(l(j))
         end do

c     Scalar case
         if (itflag.ne.2) then
            do 600 in=1,nn
               do 620 j=1,l0
                  ctnorm=dble(l(j)*(l(j)-1.0d0)*(l(j)+1)*(l(j)+2))
                  cl(j,in)=2.0d0*cl(j,in)*dble(l(j)*(l(j)+1))
                  cpl(j,in)=2.0d0*ctnorm*cpl(j,in)
     &                 *dble(l(j)*(l(j)+1))
                  ccl(j,in)=2.0d0*sqrt(ctnorm)*ccl(j,in)
     &                 *dble(l(j)*(l(j)+1))
                  ckkl(j,in)=2.0d0*ckkl(j,in)
     &                 *dble(l(j)*(l(j)+1))
                  ctkl(j,in)=2.0d0*ctkl(j,in)
     &                 *dble(l(j)*(l(j)+1))

c     TESTING. Output cls that are really calculated, from which
c     interpolation table is constructed.
c                  write(177,*)l(j),cl(j,in),cpl(j,in)

 620           continue
 600        continue

c     Making the interpolation tables to get other l-values.
            llo=1
            cllo=1.0d30
            clhi=1.0d30
            do 630 in=1,nn
               call spline(xl,cl(1,in),l0,cllo,clhi,clpr(1,in))
               call spline(xl,cpl(1,in),l0,cllo,clhi,cplpr(1,in))
               call spline(xl,ccl(1,in),l0,cllo,clhi,cclpr(1,in))
               call spline(xl,ckkl(1,in),l0,cllo,clhi,ckklpr(1,in))
               call spline(xl,ctkl(1,in),l0,cllo,clhi,ctklpr(1,in))
 630        continue
         end if
c
c
c     Tensor Case
         if (itflag.ne.0) then
c     Normalization
            do 610 j=1,l0
               ctnorm=dble(l(j)*(l(j)-1.0d0)*(l(j)+1)*(l(j)+2))
               do 640 in=1,nn
                  ctl(j,in)=ctnorm*ctl(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctel(j,in)=ctel(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctbl(j,in)=ctbl(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctcl(j,in)=sqrt(ctnorm)*ctcl(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0
 640           continue
 610        continue

c     Making the interpolation tables to get other l-values.
            cllo=1.0d30
            clhi=1.0d30
            do 650 in=1,nn
               call spline(xl,ctl(1,in),l0,cllo,clhi,ctlpr(1,in))
               call spline(xl,ctel(1,in),l0,cllo,clhi,ctelpr(1,in))
               call spline(xl,ctbl(1,in),l0,cllo,clhi,ctblpr(1,in))
               call spline(xl,ctcl(1,in),l0,cllo,clhi,ctclpr(1,in))
 650        continue
         end if

c
c     Calculating Cls for every l.
         do 700 in=1,nn
            llo=1
            do 710 il=2,l(l0)
               xi=il
               if ((xi.gt.xl(llo+1)).and.(llo.lt.l0)) then
                  llo=llo+1
               end if
               lhi=llo+1
               ho=xl(lhi)-xl(llo)
               a0=(xl(lhi)-xi)/ho
               b0=(xi-xl(llo))/ho
               if (itflag.ne.2) then
                  clint=a0*cl(llo,in)+b0*cl(lhi,in)+((a0**3-a0)
     &                 *clpr(llo,in)+(b0**3-b0)
     $                 *clpr(lhi,in))*ho*ho/6.d0
                  cplint=a0*cpl(llo,in)+b0*cpl(lhi,in)+((a0**3-a0)
     &                 *cplpr(llo,in)+(b0**3-b0)
     $                 *cplpr(lhi,in))*ho*ho/6.d0
                  cclint=a0*ccl(llo,in)+b0*ccl(lhi,in)+((a0**3-a0)
     &                 *cclpr(llo,in)+(b0**3-b0)
     $                 *cclpr(lhi,in))*ho*ho/6.d0
                  ckklint=a0*ckkl(llo,in)+b0*ckkl(lhi,in)+((a0**3-a0)
     &                 *ckklpr(llo,in)+(b0**3-b0)
     $                 *ckklpr(lhi,in))*ho*ho/6.d0
                  ctklint=a0*ctkl(llo,in)+b0*ctkl(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ctklpr(llo,in)+(b0**3-b0)
     $                 *ctklpr(lhi,in))*ho*ho/6.d0
               end if
               if (itflag.ne.0) then
                  ctlint=a0*ctl(llo,in)+b0*ctl(lhi,in)+((a0**3-a0)
     &                 *ctlpr(llo,in)+(b0**3-b0)
     $                 *ctlpr(lhi,in))*ho*ho/6.d0
                  ctelint=a0*ctel(llo,in)+b0*ctel(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ctelpr(llo,in)+(b0**3-b0)
     $                 *ctelpr(lhi,in))*ho*ho/6.d0
                  ctblint=a0*ctbl(llo,in)+b0*ctbl(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ctblpr(llo,in)+(b0**3-b0)
     $                 *ctblpr(lhi,in))*ho*ho/6.d0
                  ctclint=a0*ctcl(llo,in)+b0*ctcl(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ctclpr(llo,in)+(b0**3-b0)
     $                 *ctclpr(lhi,in))*ho*ho/6.d0
               end if

#ifdef UNNORM
               if (itflag.ne.2) then
                  clts(il,in)=clint
                  cles(il,in)=cplint
                  clcs(il,in)=cclint
                  clkk(il,in)=ckklint
                  cltk(il,in)=ctklint
               end if
               if (itflag.ne.0) then
                  cltt(il,in)=ctlint
                  clet(il,in)=ctelint
                  clbt(il,in)=ctblint
                  clct(il,in)=ctclint
               end if
#else
               if (itflag.ne.2) then
                  clts(il,in)=clint/cl(1,in)
                  cles(il,in)=cplint/cl(1,in)
                  clcs(il,in)=cclint/cl(1,in)
                  clkk(il,in)=ckklint/cl(1,in)
                  cltk(il,in)=ctklint/cl(1,in)
               end if
               if (itflag.ne.0) then
                  cltt(il,in)=ctlint/ctl(1,in)
                  clet(il,in)=ctelint/ctl(1,in)
                  clbt(il,in)=ctblint/ctl(1,in)
                  clct(il,in)=ctclint/ctl(1,in)
               end if
#endif

 710        continue
            if (itflag.ne.2) then
               clts(1,in)=cl(1,in)
            end if
            if (itflag.ne.0) then
               cltt(1,in)=ctl(1,in)
            end if
 700     continue
c
      end if

c     TIMING
c      timeprev=actual
c      actual=etime(tarray)
c      write(50,*)actual-timeprev,' END FINAL, NOW BACK'
c      write(*,*)actual-timeprev,' END FINAL, NOW BACK'

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oderivs(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
        implicit double precision (a-h,o-z)
        integer initfl
        dimension y(n),yprime(n)
        integer ndyn

c
        common /qparm/ wdyn,ndyn
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /cosmoparm/ beta,beta2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /initcase/ initfl
        common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
        parameter (lmaxnu0=25,nqmax0=15)
        parameter (lmx0=30)
        parameter (ep0=1.0d-2)
c  Internal variables.
        dimension denl(lmx0),dlfdlq(nqmax0),akv(nqmax0)
        dimension b(lmx0+1)
        common /store/ denl,dlfdlq
c
        common /epsilon/ epsw
        common /curvature/ curv,r,kcurv
        common /curvs/ aux,b
c
c ep is used to stop the tight coupling approximation.

        ak=sqrt(beta2-curv)
        ak2=ak*ak
        if (ak.gt.0.06d0*epsw) then
           ep=ep0
        else
           ep=1.1d0*ep0
        end if
        tau=x
        a=y(1)
c       ahdot=y(2)
        eta=y(3)
c  CDM.
        deltac=y(4)
        thetac=y(5)
c  Baryons.
        deltab=y(6)
        thetab=y(7)
c  Photons.
        deltag=y(8)
        thetag=y(9)
        shearg=y(10)/2.0d0
c  Polarization term.
        polter=y(10)+y(9+lmax)+y(11+lmax)
c  Massless neutrinos.
        deltar=y(10+2*lmax)
        thetar=y(11+2*lmax)
        shearr=y(12+2*lmax)/2.0d0
c
        a2=a*a
        call thermo(tau,cs2,opac,dlxedla)
c Tight Coupling parameters
        tcp=0.0d0
c       tcp1=ak/opac
        tcp1=beta/opac
        tcp2=1.0d0/(opac*tau)
c  Photon mass density over baryon mass density.
        photbar=grhog/(grhom*omegab*a)
        pb43=4.0d0/3.0d0*photbar

c Tight Coupling parameters
        tcp=0.0d0
        tcp1=ak/opac
        tcp2=1.0d0/(opac*tau)
        tcp3=ak*tau

        tcpa=0.0d0
        tcpb=0.0d0
        epstc=dmax1(tcp1,tcp2)
        if ((epstc.gt.5.0d-3).and.(a.gt.1.35d-5)) tcpa=1.0d0
        epstc=epstc/(1.0d0+pb43)
        if ((epstc.gt.5.0d-3).and.(a.gt.1.35d-4)) tcpb=1.0d0

c  Compute expansion rate.
        if (amnu.eq.0.0d0) then
           rhonu=1.0d0
           pnu=1.0d0/3.0d0
           drhonu=0.0d0
           fnu=0.0d0
           dpnu=0.0d0
           shearnu=0.0d0
        else
           call nu1(a,rhonu,pnu)
           call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
        end if

                omegavdyn = omegav * dynrho(a)
                weos = wdyn_func(a)

c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
        grho=grhom*(omegac+omegab)/a+(grhog+grhor
     $       *annur+grhonr*annunr*rhonu)/a2
     2      +grhom*omegavdyn*a2+grhom*omegak

#ifdef DIM

c energy density for  5 dimensions
        grho=(sqrt(grho-grhom*omegak)+sqrt(omegav*grhom*a2))**2
     2  +grhom*omegak

#endif

        adotoa=sqrt(grho/3.0d0)
        yprime(1)=adotoa*a
        gpres=((grhog+grhor*annur)/3.0d0+grhonr*annunr*pnu)/a2
     2      +weos*grhom*omegavdyn*a2

        if (ndyn.eq.1.or.ndyn.eq.2) then

        rc_phi = y(n-1)
        rc_psi = y(n)
         call dyn_nrg(a,grho,gpres,rc_phi,rc_psi,
     2          dgrho_phi,dgprs_phi,dgtheta_phi)



        endif

c  Evaluate metric and massive neutrino perturbations.
c       deltan=drhonu/rhonu
c       thetan=ak*fnu/(rhonu+pnu)
c  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.
        dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*annur*deltar+grhonr*annunr*drhonu)/a2

        if (ndyn.eq.1.or.ndyn.eq.2) dgrho=dgrho+dgrho_phi

        dgpres=(grhog*deltag+grhor*annur*deltar)/a2/3.0d0
     2    +grhonr*annunr*dpnu/a2

        if (ndyn.eq.1.or.ndyn.eq.2) dgpres=dgpres+dgprs_phi

c  Add a seed if desired.
        if (initfl.eq.4) dgrho=dgrho+grhom/a
        dahdotdtau=-(dgrho+3.0d0*dgpres)*a
        yprime(2)=dahdotdtau
c  Force energy conservation.
        hdot=(2.0d0*ak2*eta+dgrho)/adotoa
c  8*pi*G*(rho+P)*theta*a**2.
        dgtheta=grhom*(omegac*thetac+omegab*thetab)/a
     2    +4.0d0/3.0d0*(grhog*thetag+annur*grhor*thetar)/a2
     3    +annunr*grhonr*ak*fnu/a2

        if (ndyn.eq.1.or.ndyn.eq.2) dgtheta=dgtheta+dgtheta_phi

        etadot=0.5d0*(dgtheta+curv*hdot)/ak2
        yprime(3)=etadot

        if (tcpa.ne.1) then
c  First order approximation for photon shear
           shearg=1.0d0/opac/9.0d0*
     &      (8.0d0/3.0d0*thetag+4.0d0/3.0d0*hdot+8.0d0*etadot)
        end if

c  8*pi*G*(rho+P)*sigma*a**2.
        dgshear=4.0d0/3.0d0*(grhog*shearg+annur*grhor*shearr)/a2
     2    +annunr*grhonr*shearnu/a2
c  CDM equations of motion.
        deltacdot=-thetac-0.5d0*hdot
        yprime(4)=deltacdot
        thetacdot=-adotoa*thetac
        yprime(5)=thetacdot
c  Baryon equations of motion.
        deltabdot=-thetab-0.5d0*hdot
        yprime(6)=deltabdot
c  Need photon perturbation for first-order correction to tightly-coupled
c  baryon-photon approximation.
        deltagdot=4.0d0/3.0d0*(-thetag-0.5d0*hdot)
        drag=opac*(thetag-thetab)
        if (tcpb.eq.1) then
c  Treat baryons and photons as uncoupled.
          thetabdot=-adotoa*thetab+ak2*cs2*deltab+pb43*drag
        else
c  Treat baryons and photons as tightly coupled.
c  Zeroth-order approximation to baryon velocity.
          thetabdot=(-adotoa*thetab+ak2*cs2*deltab
     2  +ak2*pb43*(0.25d0*deltag-shearg))/(1.0d0+pb43)
c  (\ddot a)/a.

#ifdef DIM
        if (amnu.eq.0.0) then
           rhonudot=0.0d0
           shearnudot=0.0d0
        else
           call nuder(a,adotoa,rhonu,rhonudot,shearnudot,
        2                       y(iq2),yprime(iq2))
        end if
         sgrhooa2=sqrt(grho/a2)
         rgrho=sqrt(grhom*omegav+grhom*omega5
     $        +grhom*(omegab+omegac)/(a2*a)+
     $        (grhog+grhor*(annur+annunr*rhonu))/(a2*a2))
         adotdota=(sgrhooa2/rgrho*(grhom*(omegab+omegac)*
     $        (-3.0d0/(2.0d0*a))
     $        +(grhog+grhor
     $       *(annuqr+annunr*rhonu))*(-2.0d0/a2)
     $        +grhor*annunr*rhonudot/(2.0d0*a2*adotoa))+
     $        +2.0d0*grho-3.0d0*grhom*omegak)/3

#else
          adotdota=0.5d0*(adotoa*adotoa-gpres)
#endif
c  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
          slip=(2.0d0*pb43/(1.0d0+pb43)+dlxedla)
     $         *adotoa*(thetab-thetag)
     2     +1.0d0/opac*(-adotdota*thetab-adotoa*ak2*0.5d0*deltag
     3     +ak2*(cs2*deltabdot-0.25d0*deltagdot))/(1.0d0+pb43)
c  First-order approximation to baryon velocity.
          thetabdot=thetabdot+pb43/(1.0d0+pb43)*slip
        end if
        yprime(7)=thetabdot
c  Photon total intensity and polarization equations of motion.
        yprime(8)=deltagdot
        thetagdot=(-thetabdot-adotoa*thetab+ak2*cs2*deltab)/pb43
     2   +ak2*(0.25d0*deltag-b(2)*shearg)
        yprime(9)=thetagdot
        if (tcpa.eq.1) then
c  Treat baryons and photons as uncoupled.
          yprime(10)=8.0d0/15.0d0*b(2)*thetag-0.6d0*ak*b(3)*y(11)
     2        -opac*y(10)
     2        +4.0d0/15.0d0*(hdot+6.0d0*etadot*aux)*b(2)
     2        +0.1d0*opac*polter
c  Polarization equations for l = 0, 1, 2.
          yprime(9+lmax)=-ak*y(10+lmax)
     $         -opac*y(9+lmax)+0.5d0*opac*polter
          yprime(10+lmax)=ak/3.0d0*(y(9+lmax)-2.0d0*b(2)*y(11+lmax))
     2             -opac*y(10+lmax)
          yprime(11+lmax)=ak*(0.4d0*b(2)*y(10+lmax)-0.6d0*b(3)
     2             *y(12+lmax))-opac*y(11+lmax)+0.1d0*opac*polter
            do 10 l=3,lmax-1
            yprime(8+l)=ak*denl(l)*(l*b(l)
     $              *y(7+l)-(l+1)*b(l+1)*y(9+l))
     2                  -opac*y(8+l)
            yprime(9+lmax+l)=ak*denl(l)*(l*b(l)
     $           *y(8+lmax+l)-(l+1)*b(l+1)
     2                 *y(10+lmax+l))-opac*y(9+lmax+l)
10        continue
        else
c  Treat baryons and photons as tightly coupled (with no polarization).
          yprime(10)=0.0d0
          yprime(9+lmax)=0.0d0
          yprime(10+lmax)=0.0d0
          yprime(11+lmax)=0.0d0
          do 15 l=3,lmax-1
             yprime(8+l)=0.0d0
             yprime(9+lmax+l)=0.0d0
15        continue
        end if
c  Truncate moment expansion.
c       yprime(8+lmax)=ak*lmax*y(7+lmax)/(2*lmax+1)-opac*y(8+lmax)
c       yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)-opac*y(9+2*lmax)

        xc=tau/r
        alpha=ak*coshk(xc)/(sinhk(xc)*beta*r)


        if(ndyn.eq.1.or.ndyn.eq.2) then

        call dyn_phi(a,hdot,grho,gpres,rc_phi,rc_psi,
     2  rc_dphi,rc_dpsi)         
         yprime(n-1) = rc_dphi
         yprime(n)   = rc_dpsi

        endif

        yprime(8+lmax)=ak*b(lmax)*y(7+lmax)-(lmax+1)*alpha*y(8+lmax)
     2                -opac*y(8+lmax)
        yprime(9+2*lmax)=ak*b(lmax)*y(8+2*lmax)-(lmax+1)*alpha
     2                *y(9+2*lmax)-opac*y(9+2*lmax)

c  Massless neutrino equations of motion.
        deltardot=4.0d0/3.0d0*(-thetar-0.5d0*hdot)
        yprime(10+2*lmax)=deltardot
        thetardot=ak2*(0.25d0*deltar-b(2)*shearr)
        yprime(11+2*lmax)=thetardot
        yprime(12+2*lmax)=8.0d0/15.0d0*b(2)*thetar-0.6d0*ak*b(3)
     2             *y(13+2*lmax)+4.0d0/15.0d0*(hdot
     2             +6.0d0*etadot*aux)*b(2)
          do 20 l=3,lmaxnr-1
          yprime(10+2*lmax+l)=ak*denl(l)*(l*b(l)*y(9+2*lmax+l)
     2                    -(l+1)*b(l+1)*y(11+2*lmax+l))
20      continue
c  Truncate moment expansion.
c       yprime(10+2*lmax+lmaxnr)=ak*lmax
c     $       *y(9+2*lmax+lmaxnr)/(2*lmaxnr+1)

        yprime(10+2*lmax+lmaxnr)=ak*b(lmaxnr)*y(9+2*lmax+lmaxnr)
     2              -(lmaxnr+1)*alpha*y(10+2*lmax+lmaxnr)

c  Massive neutrino equations of motion.
        if (nqmax.eq.0) return
c       dq=qmax/nqmax
        dq=1.0d0
          do i=1,nqmax
          q=i*dq-0.5d0
          aq=a*amnu/q
          v=1.0d0/sqrt(1.0d0+aq*aq)
          akv(i)=ak*v
        end do
c  l = 0, 1, 2,lmaxnu.
          do 30 i=1,nqmax
          ind=iq0+i-1
          yprime(ind)=-akv(i)*y(ind+nqmax)+hdot*dlfdlq(i)/6.0d0
          ind=iq1+i-1
          yprime(ind)=akv(i)*(y(ind-nqmax)-2.0d0*b(2)
     2                                   *y(ind+nqmax))/3.0d0
          ind=iq2+i-1
          yprime(ind)=akv(i)*(2.0d0*b(2)*y(ind-nqmax)-3.0d0
     2           *b(3)*y(ind+nqmax))/5.0d0
     2           -(hdot/15.0d0+2.0d0/5.0d0*etadot*aux)
     $         *b(2)*dlfdlq(i)
          ind=10+2*lmax+lmaxnr+i+lmaxnu*nqmax
c  Truncate moment expansion.
c         yprime(ind)=akv*lmaxnu*y(ind-nqmax)/(2.0d0*lmaxnu+1)
c         yprime(ind)=akv(i)*y(ind-nqmax)-(lmaxnu+1)/tau*y(ind)
          yprime(ind)=akv(i)*y(ind-nqmax)*b(lmaxnu)
     2    -(lmaxnu+1)*y(ind)*alpha
30      continue
          do 50 l=3,lmaxnu-1
            do 40 i=1,nqmax
            ind=10+2*lmax+lmaxnr+i+l*nqmax
            yprime(ind)=akv(i)*denl(l)*(l*b(l)*y(ind-nqmax)
     &       -(l+1)*b(l+1)*y(ind+nqmax))
40        continue
50      continue
c
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ooutput(n,y,yprime,j,tau0,tau,d,dp,dk,phi)
        implicit double precision (a-h,o-z)
c
        include 'cmbfast.inc'

c       parameter (nstep0=2400)
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                                        nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
        dimension y(n),yprime(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        save /visib/

        parameter (lmx0=30)
        dimension b(lmx0+1)
        common /curvature/ curv,r,kcurv
        common /curvs/ aux,b
        common /maxreion/ taurmax,armax,jrmax

        ak=sqrt(beta2-curv)
        ak2=ak*ak

        chi=(tau0-tau)/r
        rsinh2=(r*sinhk(chi))**2
        a=y(1)
        a2=a*a
c       ahdot=y(2)
        eta=y(3)
        etadot=yprime(3)
        alpha=(hdot+6.0d0*etadot*aux)/(2.0d0*ak2)
        alphadot=-3.0d0*dgshear/(2.0d0*ak2*b(2))+eta*aux-2.0d0
     2                               *adotoa*alpha
c  Baryons.
c       deltab=y(6)
        thetab=y(7)
        thetabdot=yprime(7)
c  Photons.
        deltag=y(8)
c       thetag=y(9)
c       shearg=y(10)/2.0d0
        thetagdot=yprime(9)
        sheargdot=yprime(10)/2.0d0
c  Polarization term.
        polter=y(10)+y(9+lmaxg)+y(11+lmaxg)
c        coupl=8.0d0*(b(2)*thetag+ak2*alpha*b(2))/15.0d0-ak*0.6d0
c     2              *b(3)*(y(11)+y(10+lmaxg))-ak*
c     2              (1.0d0-0.4d0*b(2))*y(12+lmaxg)
        coupldot=8.0d0*(b(2)*thetagdot+ak2*alphadot*b(2))/15.0d0
     2      -ak*0.6d0*b(3)*(yprime(11)+yprime(10+lmaxg))
     2       +ak*(1.0d0-0.4d0*b(2))*yprime(12+lmaxg)
        polterdot=yprime(10)+yprime(9+lmaxg)+yprime(11+lmaxg)
        polterddot=coupldot-0.3d0*(dopac(j)*polter
     $       +opac(j)*polterdot)
c  Massless neutrinos.
c       deltar=y(10+2*lmaxg)
c       thetar=y(11+2*lmaxg)
c       shearr=y(12+2*lmaxg)/2.0d0
        shearrdot=yprime(12+2*lmaxg)/2.0d0
c  Second derivative of expansion rate
        if (amnu.eq.0.0) then
           rhonudot=0.0d0
           shearnudot=0.0d0
        else
           call nuder(a,adotoa,rhonu,rhonudot,shearnudot,
     2                  y(iq2),yprime(iq2))
        end if

        omegavdyn = omegav * dynrho(a)
        weos = wdyn_func(a)

        grhodot=(-grhom*(omegac+omegab)/a-2.0d0
     2      *(grhog+grhor*annur+grhonr*annunr*rhonu)/a2)*adotoa
     3      +grhonr*annunr*rhonudot/a2

        grhodot=grhodot-((1.0d0+3.0d0*weos)*grhom
     $       *omegavdyn*a2)*adotoa

#ifdef DIM

        grhodot=0.0d0
        grho=(sqrt(grhom*(omegac+omegab)/a+(grhog+grhor
        $               *(annur+annunr*rhonu))/a2
        2           +grhom*omega5*a2+grhom*omegav*a2)
        $           +sqrt(omega5*grhom)*a)**2+grhom*omegak

         sgrhooa2=sqrt(grho/a2)
        rgrho=sqrt(grhom*omega5+grhom*omegav
     $        +grhom*(omegab+omegac)/(a2*a)+
     $        (grhog+grhor*(annur+annunr*rhonu))/(a2*a2))

        grhodot=sgrhooa2/rgrho*(grhom*(omegab+omegac)*(-3.0d0/a)+
     $       (grhog+grhor*(annur+annunr*rhonu))
     $       *(-4.0d0/a2)+
     $       grhor*annunr*rhonudot/(a2*adotoa))+
     $       2.0d0*grho-2.0d0*grhom*omegak

#endif

        adotoadot=grhodot/(6.0d0*adotoa)
c  Derivative of the shear
        dgsheardot=4.0d0/3.0d0*(grhog*sheargdot+annur*grhor
     2  *shearrdot)/a2-2.0d0*adotoa
     $       *dgshear+annunr*grhonr*shearnudot/a2
c  Calculation of the sources
c        alphaddot=-3.0*dgsheardot/(2.0d0*ak2)+etadot*aux
c     2      -2.0d0*adotoadot*alpha-2.0d0*adotoa*alphadot
c
c         s1=etadot*aux+alphaddot
        alphaddot=-3.0d0*dgsheardot/(2.0d0*ak2*b(2))+etadot*aux
     2      -2.0d0*adotoadot*alpha-2.0d0*adotoa*alphadot
c
         s1=etadot*aux+alphaddot
         s2=2.0d0*alphadot
c
c     TESTING
c         s1=0.0d0

        d=expmmu(j)*s1+vis(j)*(0.25d0*deltag
     2     +s2+polter/16.0d0/b(2)+thetabdot/ak2+3.0d0/16.0d0
     3     *polterddot/ak2/b(2))+dvis(j)
     $        *(alpha+thetab/ak2+3.0d0/8.0d0
     4     *polterdot/ak2/b(2))+ddvis(j)*3.0d0/16.0d0*polter
     $        /ak2/b(2)

c     TESTING. Use only delta_gamma contribution.
c        d=vis(j)*(0.25d0*deltag+alphadot)

        if (chi.gt.0.0d0) then
c           dp=vis(j)*3.0d0/16.0d0*polter/ak2/b(2)/rsinh2
c     TESTING forget about spin 2 nature of polarization
c           dp=vis(j)*3.0d0/16.0d0*polter/b(2)

c     I will Add the /ak2/rsinh2 in the time integration,
c     this is more stable, in the suns I had problems
c     with some models, presumably due to the interpolation,
c     although it was working fine in the PC

           dp=vis(j)*3.0d0/16.0d0*polter/b(2)

        else
           dp=0.0d0
        end if
c approximate epoch of recombination
c        chir=tau0*(1-1.0/sqrt(1100.0))/r
        chir=(tau0-taurmax)/r
        if (chi.lt.chir) then
c lensing convergence, expmmu is an approximation
c one should integrate visibility function over r(chi)/r(tau)
c but the error is harmless
            phi=eta-adotoa*alpha
            glens=r*sinh(chir-chi)*r*sinh(chi)/(r*sinh(chir))
            dk=glens*ak2*phi*expmmu(j)
        else
            dk=0.0d0
        end if

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oderivst(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
        implicit double precision (a-h,o-z)
        dimension y(n),yprime(n)

        include 'cmbfast.inc'

        parameter (sqrt6=2.4494897d0)

c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                                          nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     $       grhonr
      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),itflag,
     $       lmaxt
c
        parameter (lmaxt0=10)
        parameter (ep0=1.0d-2)
c  Internal variables.
        dimension bt1(lmaxt0+1),bt2(lmaxt0+1)
c
        common /epsilon/ epsw
        common /curvature/ curv,r,kcurv
        common /curvt/ bt1,bt2

c
c ep is used to stop the tight coupling approximation.

        ak=sqrt(beta2-3.0d0*curv)
        ak2=ak*ak
c       if (ak.gt.epsw) then
        if (ak.gt.0.06d0*epsw) then
           ep=ep0
        else
           ep=1.11d0*ep0
        end if

        tau=x
        a=y(1)

        a2=a*a
        call thermo(tau,cs2,opac,dlxedla)


c Tight Coupling parameters
        tcp=0.0d0
c       tcp1=ak/opac
        tcp1=beta/opac

        tcp2=1.0d0/(opac*tau)
        if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
           tcp=1.0d0
        end if

c  Compute expansion rate.
        if (amnu.eq.0.0d0) then
           rhonu=1.0d0
           pnu=1.0d0/3.0d0
           drhonu=0.0d0
           fnu=0.0d0
           dpnu=0.0d0
           shearnu=0.0d0
        else
           call nu1(a,rhonu,pnu)
           call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
        end if

         omegavdyn = omegav * dynrho(a)

c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
        grho=grhom*(omegac+omegab)/a+(grhog+grhor
     $       *annur+grhonr*annunr*rhonu)/a2
     2      +grhom*omegavdyn*a2+grhom*omegak

#ifdef DIM
c energy density for  5 dimensions
        grho=(sqrt(grho-grhom*omegak)+sqrt(omegav*grhom*a2))**2
     2  +grhom*omegak

#endif

        adotoa=sqrt(grho/3.0d0)
        yprime(1)=adotoa*a

c Tensors
        ht=y(2)
        htpr=y(3)
        yprime(2)=htpr
c add the source
        htdpr=-2*adotoa*htpr-(ak2+2.0d0*curv)*ht
        yprime(3)=htdpr

c Photon Perturbations

c Tight Coupling parameters
        tcp=0.0d0
        tcp1=ak/opac
        tcp2=1.0d0/(opac*tau)

        if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
           tcp=1.0d0
        end if



        ind1=2
        ind2=ind1+lmaxt-1
        ind3=ind2+lmaxt-1
        p2=0.1d0*(y(ind1+2)-sqrt6*y(ind2+2))
        sourcet=opac*p2-htpr
        sourcee=-opac*sqrt6*p2
        sourceb=0.0d0


        if (tcp.eq.1) then


c No tight coupling.
c Temperature
           yprime(ind1+2)=beta*(
     2     -bt1(3)*y(ind1+3)/7.0d0)-opac*y(ind1+2)
     3      + sourcet
c E polarization
           yprime(ind2+2)=beta*(
     2     -2.0d0*y(ind3+2)/3.0d0
     3     -bt2(3)/7.0d0*
     4     y(ind2+3))-opac*y(ind2+2)
     5     + sourcee
c B polarization
           yprime(ind3+2)=beta*(
     2     +2.0d0*y(ind2+2)/3.0d0
     3     -bt2(3)/7.0d0*
     4     y(ind3+3))-opac*y(ind3+2)
     5     + sourceb


          do l=3,lmaxt-1
c Temperature
           yprime(ind1+l)=beta*(bt1(l)*y(ind1+l-1)/(2.0d0*l-1.0d0)
     2     -bt1(l+1)*y(ind1+l+1)/(2.0d0*l+3.0d0))-opac*y(ind1+l)
c E polarization
           yprime(ind2+l)=beta*(bt2(l)/(2.0d0*l-1.0d0)
     2      *y(ind2+l-1)-4.0d0/l/(l+1.0d0)*y(ind3+l)
     3     -bt2(l+1)/(2.0d0*l+3.0d0)*
     4     y(ind2+l+1))-opac*y(ind2+l)
c B polarization
           yprime(ind3+l)=beta*(bt2(l)/(2.0d0*l-1.0d0)*y(ind3+l-1)
     2     +4.0d0/l/(l+1.0d0)*y(ind2+l)
     3     -bt2(l+1)/(2.0d0*l+3.0d0)*
     4     y(ind3+l+1))-opac*y(ind3+l)
          enddo

          xc=tau/r
          betar=beta*r
          alpha=coshk(xc)/sinhk(xc)

c close the hierarchy
c Just set to zero last multipole. Fancier scheme is not working well.
c This is faster and still accurate.

c close the hierarchy
c T
        yprime(lmaxt+ind1)=0.0d0
c E
        yprime(ind2+lmaxt)=0.0d0
c B
        yprime(ind3+lmaxt)=0.0d0

        else
c tight coupling

         shearg=-4.0d0*htpr/opac/3.0d0
         e2=-sqrt(6.0d0)/4.0d0*shearg
         b2=0.0d0
         sheargdot=0.0d0
         e2dot=0.0d0
         b2dot=0.0d0
         y(ind1+2)=shearg
         y(ind2+2)=e2
         y(ind3+2)=b2
c rest is 0
         yprime(ind1+2)=sheargdot
         yprime(ind2+2)=e2dot
         yprime(ind3+2)=b2dot
         do l=3,lmaxt-1
          yprime(l+ind1)=0.0d0
          yprime(ind2+l)=0.0d0
          yprime(ind3+l)=0.0d0
         enddo
        end if

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ooutputt(n,y,yprime,j,tau0,tau,dt,dte,dtb)
        implicit double precision (a-h,o-z)
c
        include 'cmbfast.inc'

c       parameter (nstep0=2400)
        common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                                        nqmax,iq0,iq1,iq2
c
        dimension y(n),yprime(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        save /visib/

        parameter (sqrt6=2.4494897d0)
        parameter (lmaxt0=10)

        dimension bt1(lmaxt0+1),bt2(lmaxt0+1)
        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $       itflag,lmaxt
        common /curvature/ curv,r,kcurv
        common /curvt/ bt1,bt2
c

c       ak=sqrt(beta2-3.0d0*curv)

        chi=(tau0-tau)/r
        sinhchi=sinhk(chi)
        rsinh=r*sinhchi
        coshchi=coshk(chi)
        coshchi2=coshchi**2
        rsinh2=rsinh**2
c       a=y(1)

c Tensors
c
        ind1=2
        ind2=ind1+lmaxt-1
        ind3=ind2+lmaxt-1
        htpr=y(3)
        htdpr=yprime(3)

        p2=0.1d0*(y(ind1+2)-sqrt6*y(ind2+2))

        p2dot=0.1d0*(yprime(ind1+2)-sqrt6*yprime(ind2+2))

c        sourcet=-opac*p2-htpr
        sourcetdot=opac(j)*p2dot+dopac(j)*p2-htdpr

c        sourcee=sqrt6*opac(j)*p2
        sourceedot=-(sqrt6*dopac(j)*p2+sqrt6*opac(j)*p2dot)

        t2ddot=beta*(
     2     -bt1(3)*yprime(ind1+3)/7.0d0)-opac(j)*yprime(ind1+2)
     3     -dopac(j)*y(ind1+2) + sourcetdot
        e2ddot=beta*(
     2     -2.0d0*yprime(ind3+2)/3.0d0
     3     -bt2(3)/7.0d0*
     4     yprime(ind2+3))-opac(j)*yprime(ind2+2)
     5     -dopac(j)*y(ind2+2)+ sourceedot

        p2ddot=0.1d0*(t2ddot-sqrt6*e2ddot)
c

        betar=beta*r
        abeta=r*r/sqrt(betar**2-kcurv*4.0d0)
     2             /sqrt(betar**2-kcurv)

        x=beta*rsinh
        if (x.gt.1.0d-8) then

           dt=abeta*(-expmmu(j)*htpr+vis(j)*p2)/rsinh2
           dte=(ddvis(j)*p2+2.0d0*dvis(j)*p2dot+vis(j)*p2ddot)
     &          +4.0d0*coshchi*(dvis(j)*p2+vis(j)*p2dot)/rsinh
     &          + ( -beta2+3.0d0/r/r+6.0d0*coshchi2/rsinh2)
     $          *vis(j)*p2

           dte=dte*abeta


           dtb=(dvis(j)*p2+vis(j)*p2dot+2.0d0*coshchi/rsinh
     &          *vis(j)*p2)
           dtb=2.0d0*beta*abeta*dtb

        else
           dt=0.0d0
           dte=0.0d0
           dtb=0.0d0
        end if
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oinithermo(taumin,taumax,tau0,taurend,
     2      dtau0,nstep)
c  Compute and save unperturbed baryon temperature and ionization fraction
c  as a function of time.  With nthermo=10000, xe(tau) has a relative
c accuracy (numerical integration precision) better than 1.e-5.
c
        implicit double precision (a-h,o-z)
c
        include 'cmbfast.inc'

        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr

        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
c
        parameter (barssc0=9.1820d-14)
c       parameter (nstep0=2400)
c
        double precision atemp(nthermo),datemp(nthermo)
        dimension tb(nthermo),cs2(nthermo),xe(nthermo)
        dimension dcs2(nthermo)
        dimension dotmu(nthermo),sdotmu(nthermo),emmu(nthermo)
        dimension demmu(nthermo),ddotmu(nthermo)
        dimension dddotmu(nthermo),ddddotmu(nthermo)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        double precision atau0(nstep0)

        integer nreg(20),nr,rcflag
        double precision dtaureg(20)

        common /recfl/rcflag
        common /timstp1/ irec, nr,nreg,itflag
        common /timstp2/dtaureg,atau0
        common /maxreion/ taurmax,armax,jrmax
        common /thermod/ tauminn,dlntau,dotmu,ddotmu,cs2,dcs2
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
        save /thermod/
        save /visib/

        dimension dlxedla(nthermo),ddlxedla(nthermo)
        common /thermod2/ dlxedla,ddlxedla
        save /thermod2/

        common /tnowkb/twkb1,twkb2,twkb3,twkb4
c
        ncount=0
        thomc0=5.0577d-8*tcmb**4
        akthom=(2.3048d-9)*(1-yhe)*omegab*h0*h0
        tauminn=0.05d0*taumin
        dlntau=log(tau0/tauminn)/(nthermo-1)
c
c  Initial conditions: assume radiation-dominated universe.
        tau01=tauminn
        adot0=adotrad
        astart=adotrad*tauminn
        a0=astart
        a02=a0*a0
c  Assume that any entropy generation occurs before tauminn.
c  This gives wrong temperature before pair annihilation, but
c  the error is harmless.
        tb(1)=tcmb/a0
        xe0=1.0d0
        x1=0.0d0
        x2=1.0d0
        xe(1)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
        barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(1))
        cs2(1)=4.0d0/3.0d0*barssc*tb(1)
        dotmu(1)=xe(1)*akthom/a02
        sdotmu(1)=0
        atemp(1)=astart
        dlxedla(1)=0.0d0
c

          do 10 i=2,nthermo
          tau=tauminn*exp((i-1)*dlntau)
          dtau=tau-tau01
c  Integrate Friedmann equation using inverse trapezoidal rule.
          a=a0+adot0*dtau
          a2=a*a
          call nu1(a,rhonu,pnu)

         omegavdyn = omegav * dynrho(a)

          grho=grhom*(omegac+omegab)/a+(grhog+grhor
     $         *annur+grhonr*annunr*rhonu)
     2      /a2+grhom*omegavdyn*a2+grhom*omegak

#ifdef DIM
c energy density for  5 dimensions
        grho=(sqrt(grho-grhom*omegak)+sqrt(omegav*grhom*a2))**2
     2  +grhom*omegak

#endif
          adot=sqrt(grho/3.0d0)*a
          a=a0+2.0d0*dtau/(1.0d0/adot0+1.0d0/adot)
          atemp(i)=a
c  Baryon temperature evolution: adiabatic except for Thomson cooling.
c  Use  quadrature solution.
          tg0=tcmb/a0
          ahalf=0.5d0*(a0+a)
          adothalf=0.5d0*(adot0+adot)
c  fe=number of free electrons divided by total number of free baryon
c  particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
c  more accuracy is required (unlikely) then this can be iterated with
c  the solution of the ionization equation.
          fe=(1.d0-yhe)*xe(i-1)/(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i-1))
          thomc=thomc0*fe/adothalf/ahalf**3
          etc=exp(-thomc*(a-a0))
          a2t=a0*a0*(tb(i-1)-tg0)*etc-tcmb/thomc*(1.d0-etc)
          tb(i)=tcmb/a+a2t/(a*a)
c If there is re-ionization, smoothly increase xe to the
c requested value.
          if ((zri.ne.0.0d0).and.(tau.gt.(9.0d0*taurist
     &                              /10.0d0))) then
             if(ncount.eq.0) then
                ncount=i-1
             end if

c     SMOOTH REIONIZATION

             xod=150.0d0*(tau-taurist)/taurist
             if (xod.gt.100d0) then
                tgh=1
             else
                tgh=(exp(xod)-exp(-xod))/(exp(xod)+exp(-xod))
             end if
             xe(i)=(rif-xe(ncount))*(tgh+1.0d0)/2.0d0+xe(ncount)
             dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)
          else
c  Integrate ionization equation.
             if (rcflag.eq.0) then
                tbhalf=0.5d0*(tb(i-1)+tb(i))
                call ionize(tbhalf,ahalf,adothalf,dtau,xe0)
                call ionhe(tb(i),a,xe0,x1,x2)
                xe(i)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
                dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)
             else
                call recint(a,xe(i))
                dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)
             endif
          end if
c     TESTING
c         write(333,*)tau,xe(i),taurist

c  Baryon sound speed squared (over c**2).
          dtbdla=-2.0d0*tb(i)-thomc*adothalf/adot*(a*tb(i)-tcmb)
          barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i))
          cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3.0d0)

c Calculation of oppacity
          dotmu(i)=xe(i)*akthom/a2
c
 15       a0=a
          tau01=tau
          adot0=adot
 10     continue

      if ((xe(nthermo).lt.rif).and.(zri.ne.0.0d0)) then
         write(*,*)'Warning: We use a smooth function to '
         write(*,*)'approach youre specified reionization'
         write(*,*)'fraction. The redshift that is deduced from'
         write(*,*)'youre input paprameters is so low that our'
         write(*,*)'smooth function does not reach the required'
         write(*,*)'value. You should go in to subroutine oinithermo'
         write(*,*)'and play with the shape of this smooth function.'
         write(*,*)'Search for SMOOTH REIONIZATION for the'
         write(*,*)' place where the function is set.'
      end if

      do j1=1,nthermo
         sdotmu(j1)=0.0d0
      end do
      sdotmu(nthermo)=0.0d0
      do j1=nthermo-1,1,-1
         tmp1=dotmu(j1)*tauminn*exp((j1-1)*dlntau)
         tmp2=dotmu(j1+1)*tauminn*exp((j1)*dlntau)
         sdotmu(j1)=sdotmu(j1+1)+0.5d0*(tmp1+tmp2)*dlntau
         emmu(j1)=exp(-sdotmu(j1))
      end do

        iv=0
        vfi=0.0d0
c Getting the starting and finishing times for decoupling.
        if (ncount.eq.0) then
           cf1=1.0d0
           ns=nthermo
           else
              cf1=exp(sdotmu(nthermo)-sdotmu(ncount))
              ns=ncount
           end if
        do 30 j1=1,ns
           tau=tauminn*exp(dble(j1-1)*dlntau)
           vfi=vfi+emmu(j1)*dotmu(j1)*cf1*dlntau*tau
           if ((iv.eq.0).and.(vfi.gt.1.0d-12)) then
              taurst=9.0d0/10.0d0*tauminn*exp((j1-1)*dlntau)
              iv=1
           end if
           if ((iv.eq.1).and.(vfi.gt.0.999)) then
              taurend1=(tauminn*exp((j1-1)*dlntau))
              taurend1=max(taurend1,taurend1*sqrt(2500.0d0/
     &             (omegac+omegab)/h0**2))
              iv=2
           end if
 30     continue
        if (iv.ne.2) then
           taurend1=2.5d0*(tauminn*exp((ncount-1)*dlntau))
        end if

c     TESTING: Decrease time-steps during recombination
c      write(*,*)'Enter timestep accuracy factor'
c      read(*,*)accstep
c     A small timestep is needed because of the low ks, the cancellation
c     across the visibility function are not good otherwise
#ifdef HP
      accstep=6.0d0
#else
      accstep=3.0d0
#endif

c Calculating the timesteps during recombination.
           if (dtaurec.ne.0.0d0) then
              dtaurec=min(dtaurec,taurst/40.0d0/accstep)
           else
              dtaurec=taurst/40.0d0/accstep
           end if

           if (itflag.ne.0) dtaurec=dtaurec*0.5d0

c     TESTING
c           dtaurec=dtaurec*0.3d0

           dlntau0=0.005d0
           if (itflag.ne.0) dlntau0=0.025d-1
           taurend2=dtaurec/dlntau0

           taurend=max(taurend1,taurend2)

c     We do not want to switch to wkb inside recombination. Remember
c     times to check later.
           twkb1=taurst
           twkb2=taurend
           twkb3=0.0d0
           twkb4=0.0d0

           taurend=min(taurend,9.0d0*taurist/10.0d0)

           n1=int((taurend-taurst)/dtaurec)
           dtaurec=(taurend-taurst)/dble(n1)
           n1=n1+1


c Calculating the timesteps after recombination
           if (itflag.ne.0) then
              dtau0=dtaurec*3.5d0
c     We have increased the sampling during recombination
c     by a factor of 3 but we don't want to do it after.
              dtau0=dtau0*3.0
           else
              dtau0=taumax/500.0d0
           end if
           nreg(1)=1
           nreg(2)=n1
           dtaureg(1)=dtaurec
           taulast=taurend
           nr=1


c Correcting if there is reionization.



           if (zri.ne.0) then

              dt1=0.0d0
              dt2=0.0d0
              dt3=0.0d0

c     Accstep multiplies the timestep so that we do not decrease the
c     timesteps at reionization as we change that timestep at recombination.

              if(itflag.ne.2) then
                 dtauri=1.5d0*dtaurec*accstep
              else
                 dtauri=3.0d0*dtaurec*accstep
              end if

              n21=int((0.9d0*taurist-taurend)/dtau0)
              if (n21.gt.0) then
                 nr=nr+1
                 dt1=(0.9d0*taurist-taurend)/dble(n21)
                 dtaureg(nr)=dt1
                 nreg(nr+1)=nreg(nr)+n21
              end if
c             n2=n1+n21

              nr=nr+1
              nri0=50

              dtautemp=min(0.15d0*taurist,(tauristp-0.9d0*taurist))
              dt2=dtautemp/dble(nri0)
              dtaureg(nr)=dt2
              nreg(nr+1)=nreg(nr)+nri0


              n34=int((tauristp-(taurend+n21*dt1+
     &                               nri0*dt2))/dtauri)
              if (n34.gt.0) then
                 nr=nr+1
                 dt3=(tauristp-(taurend+n21*dt1+
     &                               nri0*dt2))/dble(n34)
                 dtaureg(nr)=dt3
                 nreg(nr+1)=nreg(nr)+n34
              end if

              taulast=(taurend+n21*dt1+
     &          nri0*dt2+n34*dt3)

           end if

           nlast=int((taumax-taulast)/dtau0)
           if (nlast.gt.0) then
              nr=nr+1
              dtau0=(taumax-taulast)/dble(nlast)
              dtaureg(nr)=dtau0
              nreg(nr+1)=nreg(nr)+nlast
           end if
           nstep=nreg(nr+1)


c     TESTING
c           write(*,*)'nstep .............', nstep
c           write(*,*),taurst,taurend,dtaurec,dtau0
c           write(50,*)taurend,tauristp,taurend+n21*dt1+nri0*dt2
c           write(50,*)taurist,0.15d0*taurist,(tauristp-taurist)
c           write(50,*)taumax,taulast,nlast,nr
c           write(50,*)n21,dt1,nri0,dt2,n34,dt3


           if (nstep.gt.nstep0) then
                write(*,*)'Sorry, the arrays were dimensioned'
                write(*,*)'for a max. of'
                write(*,*)nstep0, 'timesteps,'
                write(*,*)'The model you requested needs'
                write(*,*)nstep, 'Please make the arrays bigger'
                write(*,*)'by making nstep0 bigger'
              stop
           end if

        call splini
        call splder(cs2,dcs2,nthermo)
        call splder(dotmu,ddotmu,nthermo)
        call splder(ddotmu,dddotmu,nthermo)
        call splder(dddotmu,ddddotmu,nthermo)
        call splder(emmu,demmu,nthermo)
        call splder(dlxedla,ddlxedla,nthermo)
        call splder(atemp,datemp,nthermo)

c     Find the peak of the visibility function
c     directly in this arrays:
      vismax=0.0d0
      do i=1,nthermo
         tvis=dotmu(i)*emmu(i)
         if (tvis.gt.vismax) then
            vismax=tvis
            tau=tauminn*exp((i-1)*dlntau)
            taurmax=tau
            armax=atemp(i)
         endif
      end do
c     Test that this maximum is actually during recombination and
c     not during reionization
      if (armax.gt.0.003) then
         write(*,*)"The maximum of the visibility function"
         write(*,*)"occurs at z="
         write(*,*)(1/armax-1), " which looks odd"
         write(*,*)"If you are using the k splitting technique and"
         write(*,*)"you are computing a model with very high tau"
         write(*,*)" you should worry about the accuracy"
      end if
c     TESTING
c      write(*,*)taurmax,armax
c      do i=2,nthermo
c         tau=tauminn*exp((i-1)*dlntau)
c         tvis=dotmu(i)*emmu(i)
c         write(68,'(2E15.5)')tau,tvis
c      end do

c
c Saving tau values and the quantities needed to calculate
c the derivatives of the visibility function appearing in the sources.
c     TESTING
c        write(50,*)nr,nreg(1)
c        do i=1,nr
c           write(50,*)nreg(i+1),dtaureg(i)
c        end do

        call n2tau(atau0,taurst,dtaureg,nreg,nr,nstep0)
        nstep=nreg(nr+1)

c     TESTING
c        write(*,*)tau0
c        write(50,*)nr,nreg(1),nstep
c        do i=1,nr
c           write(50,*)nreg(i+1),dtaureg(i)
c        end do

c     TIMING
c        write(50,*)'timesteps'
c        do i=1,nr
c           write(50,'(1I6,2E15.5)')nreg(i)
c     $          ,dtaureg(i),atau0(nreg(i))
c        end do

c     TESTING
c        write(*,*)'nstep=',nstep
c        open(unit=71,file='timesteps_1.dat'
c     &       ,status='unknown',form='formatted')
c        rewind 71
c        do i=1,nstep
c           write(71,*)i,atau0(i)
c        end do
c        close(71)


c     Curvature radius
        hc=2.998d5/h0
        curv=-omegak/hc/hc
        r=1.0d0/sqrt(abs(curv))

        irec=1
        if (curv.gt.0.0d0) call symtau(r,atau0,dtaureg
     $       ,nreg,nr,irec,nstep0)
        nstep=nreg(nr+1)

c     TESTING
c        write(50,*)'nstep=',nstep
c        open(unit=72,file='timesteps_2.dat'
c     &       ,status='unknown',form='formatted')
c        rewind 71
c        do i=1,nstep
c           write(72,*)i,atau0(i)
c        end do
c        close(72)

        vismax=0.0d0
        do 40 j2=1,nstep
           tau=atau0(j2)
c     Cubic-spline interpolation.
           d=log(tau/tauminn)/dlntau+1.0d0
           i=int(d)
           d=d-i
           if (i.lt.nthermo) then
              opac(j2)=dotmu(i)+d*(ddotmu(i)
     $             +d*(3.0d0*(dotmu(i+1)-dotmu(i))
     2         -2.0d0*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1)
     3         +2.0d0*(dotmu(i)-dotmu(i+1)))))
              dopac(j2)=(ddotmu(i)+d*(dddotmu(i)
     $             +d*(3.0d0*(ddotmu(i+1)
     2         -ddotmu(i))-2.0d0*dddotmu(i)
     $             -dddotmu(i+1)+d*(dddotmu(i)
     3         +dddotmu(i+1)+2.0d0*(ddotmu(i)-ddotmu(i+1))))))/(tau
     4         *dlntau)
              ddopac=(dddotmu(i)+d*(ddddotmu(i)
     $             +d*(3.0d0*(dddotmu(i+1)
     2         -dddotmu(i))-2.0d0*ddddotmu(i)-ddddotmu(i+1)
     3         +d*(ddddotmu(i)+ddddotmu(i+1)+2.0d0*(dddotmu(i)
     4         -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2))
     5         /(tau*dlntau)**2
              expmmu(j2)=emmu(i)+d*(demmu(i)
     $             +d*(3.0d0*(emmu(i+1)-emmu(i))
     2         -2.0d0*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1)
     3         +2.0d0*(emmu(i)-emmu(i+1)))))

              vis(j2)=opac(j2)*expmmu(j2)
c     Find the maximum of the visibility function in the grid of
c     used to compute the sources. This is not used anywhere other
c     than in testing outputs. The values used for shifting the
c     spectra were computed using the full vectors of lenght nthermo.
              if (vis(j2).gt.vismax) then
                 vismax=vis(j2)
                 jrmax=j2
              endif

              dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
              ddvis(j2)=expmmu(j2)*(opac(j2)**3
     $             +3*opac(j2)*dopac(j2)+ddopac)
           else
              opac(j2)=dotmu(nthermo)
              dopac(j2)=ddotmu(nthermo)
              ddopac=dddotmu(nthermo)
              expmmu(j2)=emmu(nthermo)
              vis(j2)=opac(j2)*expmmu(j2)
              dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
              ddvis(j2)=expmmu(j2)*(opac(j2)**3+3*opac(j2)
     2             *dopac(j2)+ddopac)
           end if
c     TESTING
c           write(73,*)j2,atau0(j2),vis(j2)

 40     continue

c     TESTING
c        write(50,*)'jrmax=',jrmax

c     TESTING. Check visibility function.
c      write(*,*)'Maximum of vis'
c      write(*,*)taurmax,armax
c      do j2=1,nstep
c         write(60,'(1I6,4E15.5)')j2,atau0(j2),vis(j2)
c      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine oinitial(y,tau)
c  Initial conditions.
        implicit double precision (a-h,o-z)
        integer initfl
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /cosmoparm/ beta,beta2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /initcase/ initfl
c
        parameter (lmax0=8,lmaxnr0=25,lmaxnu0=25,nqmax0=15)
        parameter (lmx0=30)
        parameter (csp=2.998d5)
        dimension b(lmx0+1)
        common /qparm/ wdyn,ndyn
        parameter(nvar0=9+2*(lmax0+1)+(lmaxnr0+1)+nqmax0
     $       *(lmaxnu0+1))
        dimension y(nvar0)
        common /curvature/ curv,r,kcurv
        common /curvs/ aux,b
c


        curv=-omegak*h0*h0/csp/csp
        ak=sqrt(beta2-curv)
        ak2=ak*ak

        aux=1.0d0/(1.0d0-3.0d0*curv/ak2)

        do l=1,lmx0+1
           b(l)=1.0d0-curv*(l*l-1.0d0)/ak2
           if (b(l).lt.0.0) b(l)=0.0d0
           b(l)=sqrt(b(l))
        end do

        a=tau*adotrad
        a2=a*a
        call nu1(a,rhonu,pnu)
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.

         omegavdyn = omegav * dynrho(a)
         weos = wdyn_func(a)

        grho=grhom*(omegac+omegab)/a+(grhog+grhor
     $       *annur+grhonr*annunr*rhonu)/a2
     2      +grhom*omegavdyn*a2+grhom*omegak

#ifdef DIM
c energy density for  5 dimensions
        grho=(sqrt(grho-grhom*omegak)+sqrt(omegav*grhom*a2))**2
     2  +grhom*omegak

#endif

c       adotoa=sqrt(grho/3.0d0)
        gpres=((grhog+grhor*annur)/3.0d0+grhonr*annunr*pnu)/a2
     2  +weos*grhom*omegavdyn*a2
        s=grho+gpres
        fracnu=4.0d0/3.0d0*(grhor*annur+grhonr*annunr)/a2/s
c  Use yrad=rho_matter/rho_rad to correct initial conditions for
c  matter+radiation.
        yrad=grhom*(omegac+omegab)*a/
     $       (grhog+grhor*annur+grhonr*annunr*rhonu)
c
c  Choose one of the following four cases for initial conditions, or
c  add your own.  Comment out the other cases.
c
        if (initfl.eq.1) then
c-------------------------------------------------------------------------------
c  First case.
c  Isentropic ("adiabatic") initial conditions.
c normalize to zeta=1 in radiation era
#ifdef OLDNORM
        psi=-1.0d0
#else
        psi=-1.0d0/(1.5d0+2.0d0*fracnu/5.0d0)
#endif
        C=(15.0d0+4.0d0*fracnu)/20.0d0*psi
        akt2=(ak*tau)**2
        h=C*akt2*(1.0d0-0.2d0*yrad)
        eta=2.0d0*C-(5.0d0+4.0d0*fracnu)/6.0d0
     $       /(15.0d0+4.0d0*fracnu)*
     2              C*akt2*(1.0d0-yrad/3.0d0)
        eta=(eta+0.5d0*curv*h/ak2)
        f1=(23.0d0+4.0d0*fracnu)/(15.0d0+4.0d0*fracnu)
        deltac=-0.5d0*h
        deltag=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0)
        deltab=0.75d0*deltag
        deltar=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0*f1)
        thetac=0.0d0
        thetag=-C/18.0d0*akt2*akt2/tau
        thetab=thetag
        thetar=f1*thetag
        shearr=4.0d0/15.0d0*ak2/s*psi*(1.0d0+7.0d0/36.0d0*yrad)/b(2)
        ahdot=2.0d0*C*ak2*tau*a*(1.0d0-0.3d0*yrad)
c-------------------------------------------------------------------------------
        else if (initfl.eq.2) then
c  Second case.
c  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.
        delta0=1.0d0
        h=delta0*yrad*(1.0d0/(1.0d0+omegab/omegac)-0.5d0*yrad)
        deltac=delta0-0.5d0*h
c  Compensate perturbation with everything else.
        deltag=-2.0d0/3.0d0*h
        deltab=0.75d0*deltag
        deltar=deltag
        thetac=0.0d0
        thetag=-h/12.0d0*ak2*tau
        thetab=thetag
        thetar=thetag
        shearr=0.0d0
        ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
        eta=-h/6.0d0
        eta=(eta+0.5d0*curv*h/ak2)
c-------------------------------------------------------------------------------
        else if (initfl.eq.3) then
c  Third case.
c  Isocurvature baryon initial conditions: perturb only baryons as a->0.
        delta0=1.0d0
        h=delta0*yrad*(1.0d0/(1.0d0+omegac/omegab)-0.5d0*yrad)
        deltab=delta0-0.5d0*h
c  Compensate perturbation with everything else.
        deltac=-0.5d0*h
        deltag=-2.0d0/3.0d0*h
        deltar=deltag
        thetac=0.0d0
        thetag=-h/12.0d0*ak2*tau
        thetab=thetag
        thetar=thetag
        shearr=0.0d0
        ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
        eta=-h/6.0d0
        eta=(eta+0.5d0*curv*h/ak2)
c-------------------------------------------------------------------------------
        else if (initfl.eq.4) then
c  Fourth case.
c  Isocurvature seed initial conditions:everything is unperturned as a->0
        delta0=1.0d0
        h=delta0*yrad*(1.0d0/(1.0d0+omegac/omegab)-0.5d0*yrad)
c  Compensate perturbation with everything else.
        deltab=-0.5d0*h
        deltac=-0.5d0*h
        deltag=-2.0d0/3.0d0*h
        deltar=deltag
        thetac=0.0d0
        thetag=-h/12.0d0*ak2*tau
        thetab=thetag
        thetar=thetag
        shearr=0.0d0
        ahdot=adotrad*h*(1.0d0-0.5d0*yrad)
        eta=-h/6.0d0
        eta=(eta+0.5d0*curv*h/ak2)
c-------------------------------------------------------------------------------
        else
          write(*,*) 'initfl must equal 1-4! initfl=',initfl
          stop
        end if
c
        deltan=deltar
        thetan=thetar
c
        y(1)=a
        y(2)=ahdot
        y(3)=eta
c  CDM.
        y(4)=deltac
        y(5)=thetac
c  Baryons.
        y(6)=deltab
        y(7)=thetab
c  Photons (total intensity and polarization).
        y(8)=deltag
        y(9)=thetag
c       shearg=0.0d0
        y(9+lmax)=0.0d0
        y(10+lmax)=0.0d0
          do 10 l=2,lmax
          y(8+l)=0.0d0
          y(9+lmax+l)=0.0d0
10      continue
c  Massless neutrinos.
        y(10+2*lmax)=deltar
        y(11+2*lmax)=thetar
        y(12+2*lmax)=shearr*2.0d0
          do 20 l=3,lmaxnr
          y(10+2*lmax+l)=0.0d0
20      continue
c  Massive neutrinos.
        if (nqmax.eq.0) go to 50
c       dq=qmax/nqmax
        dq=1.0d0
          do 40 i=1,nqmax
          q=i*dq-0.5d0
          aq=a*amnu/q
          v=1.0d0/sqrt(1.0d0+aq*aq)
          akv=ak*v
          expq=exp(-q)
          dlfdlq=-q/(1.0d0+expq)
          y(iq0+i-1)=-0.25d0*dlfdlq*deltan
c  Divide by v to get first-order correction for neutrino mass.
          y(iq1+i-1)=-dlfdlq*thetan/akv/3.0d0
          y(iq2+i-1)=-0.5d0*dlfdlq*shearr
            do 30 l=3,lmaxnu
            ind=10+2*lmax+lmaxnr+i+l*nqmax
            y(ind)=0.0d0
30        continue
40      continue
c  Check energy constraint equation.
50      call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
        deltan=drhonu/rhonu
        thetan=ak*fnu/(rhonu+pnu)
c       shearn=shearnu/(rhonu+pnu)

CQ      quintessence: fluctuating scalar field
CQ      set initial conditions for fluctuations
CQ      in the scalar field, \delta\phi and \delta\phi'
CQ
CQ      For most cases of interest, the fluctuations are
CQ      relativistic at early times, and the perturbations
CQ      damp out. Hence, the late time perturbations
CQ      induced in the scalar field _in response_ to the
CQ      perturbations in the CDM are more important than
CQ      the initial perturbations. For simplicity, we
CQ      set the initial perturbations to zero.
CQ

       if (ndyn.eq.1.or.ndyn.eq.2) then
         rc_phi = 0.0
         rc_psi = 0.0
         y(nvar0-1) = rc_phi
         y(nvar0)   = rc_psi

         call dyn_nrg(a,grho,gpres,rc_phi,rc_psi,
     2   dgrho_phi,dgprs_phi,dgtheta_phi)
         dgrho=dgrho+dgrho_phi
       endif


        dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*annur*deltar+grhonr*annunr*drhonu)/a2
c  Add a seed if desired.
        if (initfl.eq.4) dgrho=dgrho+grhom/a
c       econ=(adotoa*ahdot/a-2.0d0*ak2*eta-dgrho)/grho
        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine oinitialt(y,tau)
c  Initial conditions.
        implicit double precision (a-h,o-z)
c
        include 'cmbfast.inc'
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,
     $       dtaurec,grhonr
      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $       itflag,lmaxt
c
        common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                   nqmax,iq0,iq1,iq2

        parameter (lmaxt0=10)
        parameter (nvar0t=3*lmaxt0)
        parameter (csp=2.998d5)
c
        dimension y(nvar0t)
        dimension bt1(lmaxt0+1), bt2(lmaxt0+1)

        common /curvature/ curv,r,kcurv
        common /curvt/ bt1, bt2

        curv=-omegak*h0*h0/csp/csp
c       ak=sqrt(beta2-3.0d0*curv)
c       ak2=ak*ak
        betar=beta*r
        do l=2,lmaxt0+1
           bt1(l)=(1.0d0/dble(l)/dble(l)+1.0d0/betar/betar)*
     &        (dble(l)*dble(l)-4.0d0)
           bt2(l)=(1.0d0/dble(l)/dble(l)+1.0d0/betar/betar)
           if (bt1(l).lt.0.0) bt1(l)=0.0d0
           if (bt2(l).lt.0.0) bt2(l)=0.0d0
           bt1(l)=sqrt(bt1(l))*dble(l)
           bt2(l)=sqrt(bt2(l))*(dble(l)*dble(l)-4.0d0)
        end do

c
        a=tau*adotrad
c
        y(1)=a
c Tensor modes
c       ht=1.0d0
c       htpr=0.0d0
        y(2)=1.0d0
        y(3)=0.0d0
        ind1=2
        ind2=ind1+lmaxt-1
        ind3=ind2+lmaxt-1
        do l=2,lmaxt
           y(ind1+l)=0.0d0
           y(ind2+l)=0.0d0
           y(ind3+l)=0.0d0
        end do
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine initujl(ll,ak,chi0,tau,y1,y2,r,tau0)
c     The bessel function is different from zero starting
c     at chi0, but we have to start the integration at recombination,
c     given by tau. So this routine advances from chi0,y1,y2 to
c     tau and also updates the values of y1 and y2.
c     If we eventually stop integrating the differential
c     equation for the ujls and just call a subroutine of the form
c     ujl(chi,beta,l,K) the need to call initujl will
c     disapear.
      implicit double precision(a-h,o-z)

c     TESTING: have a flag that sets if that particular ak is outputed
c      common /akoutput/ koutflag,kfile

      include 'cmbfast.inc'

      beta=ak*r
      beta2=beta**2
      betam1=1.0d0/beta



c Integrating the diff. equation to get to tau

      chi1=(tau0-tau)/r

c     Finding the step in the for the ujl integration.

      if (chi1.gt.chi0) then

c         delchi=0.5d0*betam1
c         delchi=0.25d0*betam1

         if ((ll.gt.1).and.(ll.lt.100)) then
            delchi=0.15d0*betam1
         else
            if (ll.lt.900) then
               delchi=0.35d0*betam1
            else
               delchi=0.5d0*betam1
            end if
         end if

         sh=sinhk(chi0)
         ch=coshk(chi0)
         ap1=dble(ll*(ll+1))

         delchi=min(delchi,0.3d0*sh/ch)


         nstep=int((chi1-chi0)/delchi)+1
         delchi=(chi1-chi0)/nstep
         chi=chi0

         do i=1,nstep

c  One step in the ujl integration

            call evalujl(chi,y1,y2,delchi,beta2,ap1)
c     TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
c            if ((ll.eq.2).and.(koutflag.eq.1)) then
c               ujl=y1/sinhk(chi)
c               write(kfile,'(5E15.5,1I6)')chi,ujl,0.0d0,0.0d0,0.0d0,1
c            end if

         end do
      end if
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine intopen1(atau0,nstart,nend,dtau,ak,l,tau0
     &  ,r,s2,sp2,sk2,ds2,dsp2,dsk2,y1,y2,out1,out2,out3)

c This subroutine integrates the source*ujl.
c It calculates ujl by integrating a second order
c differential equation from initial values for calculating ujl.
c atau0 is the array with the time where the sources are stored.
c nstart and nend are the starting and finishing values of the
c integration.
c dtau is the spacing of the timestaps (they must be equally spaced)
c s2, sp2 and sk2 are the sources. ds2, dp2, dsk2 are their second
c derivatives for the spline interpolation.
c ak is the wavevector, l is the multipole.
c tau0 is the time today.
c r is the curvature radius.

      implicit double precision(a-h,o-z)
      double precision atau0(*),s2(*),sp2(*),sk2(*)
      double precision ds2(*),dsp2(*),dsk2(*)
      integer nstart,nend
c dxsource is the spacing in x of the sampling of the source.
c dxmax is the maximun step for which the ujl integration converges
c correctly.

c     TESTING. To output sources and ujl at LSS need to know when
c     this happens.
c      common /maxreion/ taurmax,armax,jrmax
c     TESTING: have a flag that sets if that particular ak is outputed
c      common /akoutput/ koutflag,kfile


      if (nstart.eq.nend) then
         out1=0.0d0
         out2=0.0d0
         out3=0.0d0
         return
      end if

      aux1=1.0d0*r/dtau
      aux2=(tau0-atau0(nstart))/dtau + nstart

      beta=ak*r
      betam1=1.0d0/beta
      beta2=beta*beta
      betam2=betam1*betam1
      ap1=dble(l*(l+1))

      chi=(tau0-atau0(nend))/r
      sh=sinhk(chi)
      shm1=1.0d0/sh
      ujl=y1*shm1

      if ((l.gt.1).and.(l.lt.100)) then
         dchimax=0.15d0*betam1
      else
         if (l.lt.900) then
            dchimax=0.35d0*betam1
         else
            dchimax=0.5d0*betam1
         end if
      end if
      dchimax=min(dchimax,abs(0.3d0*sh/(coshk(chi)+0.1d0)))
      dchisource=dtau/r

c     TESTING
c      dchimax=dchimax/2.0

      out1=0.5d0*ujl*s2(nend)
      out2=0.5d0*ujl*sp2(nend)*betam2*shm1*shm1
      out3=0.5d0*ujl*sk2(nend)

c     TESTING
c      if (koutflag.eq.1) then
c         write(419,*)l,nend,nstart,3
c         write(419,*)dchisource,dchimax
c         write(419,*)aux1,aux2,tau0
c      end if

      if (dchisource.gt.dchimax) then

c     TESTING
c      if (koutflag.eq.1) then
c         write(419,*)'a'
c      end if

c Interpolate the source.
         delchi=dchimax
         Deltachi=(atau0(nend)-atau0(nstart))/r
         nstep=int(Deltachi/delchi)+1
         delchi=Deltachi/dble(nstep)
         dtau2=r*r*delchi*delchi

         do i=1,nstep

c One step in the ujl integration

            call evalujl(chi,y1,y2,delchi,beta2,ap1)
            shm1=1.0d0/sinhk(chi)
            ujl=y1*shm1

c Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
            b=taui-dble(is)
            a=1.0d0-b
            st=a*s2(is)+b*s2(is+1)+((a**3-a)*ds2(is)+
     &     (b**3-b)*ds2(is+1))*dtau2/6.d0
            sp=a*sp2(is)+b*sp2(is+1)+((a**3-a)*dsp2(is)+
     &     (b**3-b)*dsp2(is+1))*dtau2/6.d0
            sk=a*sk2(is)+b*sk2(is+1)+((a**3-a)*dsk2(is)+
     &     (b**3-b)*dsk2(is+1))*dtau2/6.d0

c     TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
c            if ((l.eq.2).and.(koutflag.eq.1)) then
c               write(kfile,'(5E15.5,1I6)')chi,ujl,st,sp,real(1),is
c               write(kfile+10,'(5E15.5,1I6)')chi,st
c     $              ,(aux2-is)/aux1,s2(is),taui*dtau,is
c            end if

c     TESTING, Output ujl and sources at recombination for each
c     for a particular l.
c            if ((is.eq.jrmax)
c     $           .and.(l.eq.950))
c     $           write(101,'(7E15.5)')ak,ujl
c     $           ,st,sp,beta,chi,real(1.0)

            out1=out1+st*ujl
            out2=out2+sp*ujl*betam2*shm1*shm1
            out3=out3+sk*ujl
         end do

         out1=out1-0.5d0*st*ujl
         out2=out2-0.5d0*sp*ujl*betam2*shm1*shm1
         out3=out3-0.5d0*sk*ujl
         out1=out1*delchi*r
         out2=out2*delchi*r
         out3=out3*delchi*r

      else

c     TESTING
c         if (koutflag.eq.1) then
c            write(419,*)'b',nend-1,nstart
c         end if

c Compute jl at values where the source is stored.
         delchi=dchisource

         do i=(nend-1),nstart,-1

c One step in the ujl integration

            call evalujl(chi,y1,y2,delchi,beta2,ap1)
            shm1=1.0d0/sinhk(chi)
            ujl=y1*shm1
c     TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
c            if ((l.eq.2).and.(koutflag.eq.1)) then
c               ujl=y1/sinhk(chi)
c               st=s2(i)
c               write(kfile,'(5E15.5,1I6)')chi,ujl,st,sp,real(3),i
c            end if

c     TESTING, Output ujl and sources at recombination for each
c     for a particular l.
c            st=s2(i)
c            sp=sp2(i)
c            if ((i.eq.jrmax)
c     $           .and.(l.eq.1050))
c     $           write(101,'(7E15.5)')ak,ujl
c     $           ,st,sp,beta,chi,real(2.0)

            out1=out1+s2(i)*ujl
            out2=out2+sp2(i)*ujl*betam2*shm1*shm1
            out3=out3+sk2(i)*ujl
         end do
         out1=out1-0.50d0*s2(nstart)*ujl
         out2=out2-0.50d0*sp2(nstart)*ujl*betam2*shm1*shm1
         out3=out3-0.50d0*sk2(nstart)*ujl
         out1=out1*dtau
         out2=out2*dtau
         out3=out3*dtau
      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine intopen2(atau0,nstart,nend,dtau,ak,l,tau0
     &,r,s2,sp2,sk2,ds2,dsp2,dsk2,y1,y2,out1,out2,out3)

c This subroutine integrates the source*ujl.
c It approximates ujl by an asintotic formula using y1 and y2
c to match continuously with the exact ujl.
c atau0 is the array with the time where the sources are stored.
c nstart and nend are the starting and finishing values of the
c integration.
c dtau is the spacing of the timestaps (they must be equally spaced)
c s2 and sp2 are the sources. ds2 dp2 are their second derivatives for
c the spline interpolation.
c ak is the wavevector, l is the multipole.
c tau0 is the time today.
c r is the curvature radius.

      implicit double precision(a-h,o-z)
      double precision atau0(*),s2(*),ds2(*)
      double precision sp2(*),dsp2(*)
      double precision sk2(*),dsk2(*)
      integer nstart,nend
      parameter (pi=3.14159265359d0)

      common /flip/ichiflag

c     TESTING. To output sources and ujl at LSS need to know when
c     this happens.
c      common /maxreion/ taurmax,armax,jrmax
c     TESTING: have a flag that sets if that particular ak is outputed
c      common /akoutput/ koutflag,kfile

      if (nstart.eq.nend) then
         out1=0.0d0
         out2=0.0d0
         out3=0.0d0
         return
      end if

      aux1=1.0d0*r/dtau
      aux2=(tau0-atau0(nstart))/dtau + nstart
      dtau2=dtau*dtau

      betam1=1.0d0/ak/r
      betam2=betam1*betam1
      g2=dble(l*(l+1))/(ak*r)**2

      chi=(tau0-atau0(nend))/r
      shm1=1.0d0/sinhk(chi)
      amp=betam1/(1.0d0-g2*shm1*shm1)**0.25d0

      aux01=y1/amp
      if (abs(aux01).ge.1.) then
         y1=y1*abs(amp/y1)
         aux01=y1/abs(y1)
      end if

      if (y2.ge.0d0) then
         delta=asin(aux01)
      else
         delta=pi-asin(aux01)
      end if

      dchimax=0.5d0*betam1

      dchisource=dtau/r

c     There is no need to integrate modes that oscillate
c     too much in the time the source remains constant, so return.

      if (dchimax.lt.(dchisource*0.5d0)) then
         out1=0.0d0
         out2=0.0d0
         out3=0.0d0
         return
      end if

      delchi=min(dchimax,dchisource)
      dtaup=r*delchi

      nstp=int((atau0(nend)-atau0(nstart))/dtaup)+1
      dtaup=(atau0(nend)-atau0(nstart))/dble(nstp)
      delchi=dtaup/r

      ujl=y1*shm1
      out1=0.5d0*s2(nend)*ujl
      out2=0.5d0*sp2(nend)*ujl*betam2*shm1*shm1
      out3=0.5d0*sk2(nend)*ujl

c     TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
c            beta=ak*r
c            st=s2(nend)
c            if ((l.eq.2).and.(koutflag.eq.1)) then
c
c               write(kfile,'(5E15.5,1I6)')chi,ujl,st,sp,amp*shm1
c     $              *sin(delta),nend
c            end if

      delphi=ak*dtaup

c     In models where the source has been flipped because
c     chi to recombination is more than pi/2 make sure ujl
c     has the correct parity. This means that when we switch
c     to WKB the ujl is discontineous. This does not cause significant
c     error. Instead if the parity is not correct and pi/2 falls inside
c     recombination we run into big trouble because the cancelation of the
c     velocity term does not occur (because of the integration by parts) so
c     modes with a significant velocity contribution are wrong. So we need
c     to make sure that the parity is correct. In models where pi/2 is not
c     inside recombiantion, all is irrelevant.

      if (ichiflag.eq.1) then
         beta=ak*r
         delta=beta*chi-dble(l)*pi/2.0d0
         delta=mod(delta,2.0d0*pi)
      end if
      phi=delta


      do i=1,nstp
         chi=chi+delchi
         shm1=1.0d0/sinhk(chi)
         phi=phi+delphi
         ujl=betam1*shm1*sin(phi)/(1.0d0-g2*shm1*shm1)**0.25d0


c Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
            b=taui-dble(is)
            a=1.0d0-b
            st=a*s2(is)+b*s2(is+1)+((a**3-a)*ds2(is)+
     &     (b**3-b)*ds2(is+1))*dtau2/6.d0
            sp=a*sp2(is)+b*sp2(is+1)+((a**3-a)*dsp2(is)+
     &     (b**3-b)*dsp2(is+1))*dtau2/6.d0
            sk=a*sk2(is)+b*sk2(is+1)+((a**3-a)*dsk2(is)+
     &     (b**3-b)*dsk2(is+1))*dtau2/6.d0

c     TESTING: OUTPUT TIME INTEGRAL FOR ONE l and several betas.
c            beta=ak*r
c            if ((l.eq.2).and.(koutflag.eq.1)) then
c               write(kfile,'(5E15.5,1I6)')chi,ujl,st
c     $              ,sp,real(4)
c     $              ,i
c               write(kfile+10,'(5E15.5,1I6)')chi,st
c     $              ,(aux2-is)/aux1,s2(is),taui*dtau,is
c            end if
c     TESTING, Output ujl and sources at recombination for each
c     for a particular l.
c            beta=ak*r
c            if ((is.eq.jrmax)
c     $           .and.(l.eq.950))
c     $           write(101,'(7E15.5)')ak,ujl
c     $           ,st,sp,beta,chi,real(3.0)

            out1=out1+ujl*st
            out2=out2+ujl*sp*betam2*shm1*shm1
            out3=out2+ujl*sk

      end do

      out1=out1-0.5d0*ujl*s2(nstart)
      out1=out1*dtaup
      out2=out2-0.5d0*ujl*sp2(nstart)*betam2*shm1*shm1
      out2=out2*dtaup
      out3=out3-0.5d0*ujl*sk2(nstart)
      out3=out3*dtaup

      y1=ujl*sinhk(chi)
      y2=cos(phi)

      iphase=1

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine intopen1t(atau0,nstart,nend,dtau,ak,l,tau0
     &                   ,r,s2,se2,sb2,ds2,dse2,dsb2,y1,y2
     &                   ,out1,out2,out3)

c This subroutine integrates the source*ujl.
c It calculates ujl by integrating a second order
c differential equation from initial values for calculating ujl.
c atau0 is the array with the time where the sources are stored.
c nstart and nend are the starting and finishing values of the
c integration.
c dtau is the spacing of the timestaps (they must be equally spaced)
c s2, se2 and sb2 are the sources. ds2 dp2 are their second derivatives for
c the spline interpolation.
c ak is the wavevector, l is the multipole.
c tau0 is the time today.
c r is the curvature radius.

      implicit double precision(a-h,o-z)
      double precision atau0(*),s2(*),se2(*),sb2(*)
      double precision ds2(*),dse2(*),dsb2(*)
      integer nstart,nend

c dxsource is the spacing in x of the sampling of the source.
c dxmax is the maximun step for which the ujl integration converges
c correctly.

      if (nstart.eq.nend) then
         out1=0.0d0
         out2=0.0d0
         out3=0.0d0
         return
      end if

      aux1=1.0d0*r/dtau
      aux2=(tau0-atau0(nstart))/dtau + nstart

      beta=ak*r
      betam1=1.0d0/beta
      beta2=beta*beta
      ap1=dble(l*(l+1))


      chi=(tau0-atau0(nend))/r
      sh=sinhk(chi)
      ujl=y1/sh

      if ((l.gt.1).and.(l.lt.300)) then
         if (l.lt.60) then
            dchimax=0.15d0*betam1
         else
            dchimax=0.06d0*betam1
            if (l.ge.600) dchimax=0.03d0*betam1
         end if
      else
         if (l.lt.900) then
            dchimax=0.35d0*betam1
         else
            dchimax=0.5d0*betam1
         end if
      end if

c     TESTING
c      dchimax=dchimax*0.5d0


      dchimax=min(dchimax,abs(0.3d0*sh/(coshk(chi)+0.1d0)))
      dchisource=dtau/r

      out1=0.5d0*ujl*s2(nend)
      out2=0.5d0*ujl*se2(nend)
      out3=0.5d0*ujl*sb2(nend)

      if (dchisource.gt.dchimax) then
c Interpolate the source.
         delchi=dchimax
         Deltachi=(atau0(nend)-atau0(nstart))/r
         nstep=int(Deltachi/delchi)+1
         delchi=Deltachi/dble(nstep)
         dtau2=r*r*delchi*delchi

         do i=1,nstep

c One step in the ujl integration

            call evalujl(chi,y1,y2,delchi,beta2,ap1)
            ujl=y1/sinhk(chi)


c Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
            b=taui-dble(is)
            a=1.0d0-b
            st=a*s2(is)+b*s2(is+1)+((a**3-a)*ds2(is)+
     &     (b**3-b)*ds2(is+1))*dtau2/6.d0
            se=a*se2(is)+b*se2(is+1)+((a**3-a)*dse2(is)+
     &     (b**3-b)*dse2(is+1))*dtau2/6.d0
            sb=a*sb2(is)+b*sb2(is+1)+((a**3-a)*dsb2(is)+
     &     (b**3-b)*dsb2(is+1))*dtau2/6.d0


c Add.
            out1=out1+st*ujl
            out2=out2+se*ujl
            out3=out3+sb*ujl

         end do

         out1=out1-0.5d0*st*ujl
         out2=out2-0.5d0*se*ujl
         out3=out3-0.5d0*sb*ujl
         out1=out1*delchi*r
         out2=out2*delchi*r
         out3=out3*delchi*r

      else

c Compute jl at values where the source is stored.
         delchi=dchisource
         do i=(nend-1),nstart,-1

c One step in the ujl integration

            call evalujl(chi,y1,y2,delchi,beta2,ap1)
            ujl=y1/sinhk(chi)
            out1=out1+s2(i)*ujl
            out2=out2+se2(i)*ujl
            out3=out3+sb2(i)*ujl
         end do
         out1=out1-0.50d0*s2(nstart)*ujl
         out2=out2-0.50d0*se2(nstart)*ujl
         out3=out3-0.50d0*sb2(nstart)*ujl
         out1=out1*dtau
         out2=out2*dtau
         out3=out3*dtau
      end if


      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine intopen2t(atau0,nstart,nend,dtau,ak,l,tau0
     &                   ,r,s2,se2,sb2,ds2,dse2,dsb2,y1,y2
     &                   ,out1,out2,out3)

c This subroutine integrates the source*ujl.
c It approximates ujl by an asintotic formula using y1 and y2
c to match continuously with the exact ujl.
c atau0 is the array with the time where the sources are stored.
c nstart and nend are the starting and finishing values of the
c integration.
c dtau is the spacing of the timestaps (they must be equally spaced)
c s2, se2 and sb2 are the sources. ds2 dp2 are their second derivatives for
c the spline interpolation.
c ak is the wavevector, l is the multipole.
c tau0 is the time today.
c r is the curvature radius.

      implicit double precision(a-h,o-z)
      double precision atau0(*),s2(*),ds2(*)
      double precision se2(*),dse2(*)
      double precision sb2(*),dsb2(*)
      integer nstart,nend
      parameter (pi=3.14159265359d0)

      common /flip/ichiflag

      if (nstart.eq.nend) then
         out1=0.0d0
         out2=0.0d0
         out3=0.0d0
         return
      end if

      aux1=1.0d0*r/dtau
      aux2=(tau0-atau0(nstart))/dtau + nstart
      dtau2=dtau*dtau

      betam1=1.0d0/ak/r
      g2=dble(l*(l+1))/(ak*r)**2

      chi=(tau0-atau0(nend))/r
      shm1=1.0d0/sinhk(chi)
      amp=betam1/(1.0d0-g2*shm1*shm1)**0.25d0

      aux01=y1/amp
      if (abs(aux01).ge.1.) then
         y1=y1*abs(amp/y1)
         aux01=y1/abs(y1)
      end if

      if (y2.ge.0d0) then
         delta=asin(aux01)
      else
         delta=pi-asin(aux01)
      end if

      dchimax=0.5d0*betam1

      dchisource=dtau/r

      delchi=min(dchimax,dchisource)
      dtaup=r*delchi

      nstp=int((atau0(nend)-atau0(nstart))/dtaup)+1
      dtaup=(atau0(nend)-atau0(nstart))/dble(nstp)
      delchi=dtaup/r

      ujl=y1*shm1
      out1=0.5d0*s2(nend)*ujl
      out2=0.5d0*se2(nend)*ujl
      out3=0.5d0*sb2(nend)*ujl

      delphi=ak*dtaup

      if (ichiflag.eq.1) then
         beta=ak*r
         delta=beta*chi-dble(l)*pi/2.0d0
         delta=mod(delta,2.0d0*pi)
      end if
      phi=delta

      do i=1,nstp
         chi=chi+delchi
         shm1=1.0d0/sinhk(chi)
         phi=phi+delphi
         ujl=betam1*shm1*sin(phi)/(1.0d0-g2*shm1*shm1)**0.25d0

c Interpolate the source
            taui=aux2-aux1*chi
            is=int(taui)
            b=taui-dble(is)
            a=1.0d0-b
            st=a*s2(is)+b*s2(is+1)+((a**3-a)*ds2(is)+
     &     (b**3-b)*ds2(is+1))*dtau2/6.d0
            se=a*se2(is)+b*se2(is+1)+((a**3-a)*dse2(is)+
     &     (b**3-b)*dse2(is+1))*dtau2/6.d0
            sb=a*sb2(is)+b*sb2(is+1)+((a**3-a)*dsb2(is)+
     &     (b**3-b)*dsb2(is+1))*dtau2/6.d0
            out1=out1+ujl*st
            out2=out2+ujl*se
            out3=out3+ujl*sb
      end do

      out1=out1-0.5d0*ujl*s2(nstart)
      out1=out1*dtaup
      out2=out2-0.5d0*ujl*se2(nstart)
      out2=out2*dtaup
      out3=out3-0.5d0*ujl*sb2(nstart)
      out3=out3*dtaup

      y1=ujl*sinhk(chi)
      y2=cos(phi)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine initujl0(filename)

c     This subroutine reads the ujl files from disk and
c     initializes other variables needed in CMBFAST.

      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      parameter (nbmax=1000)

      parameter (d0hi=1.0d40,d0lo=1.0d40)

      double precision abeta(nbmax),ax(nbmax,lmax)
      double precision ay2(nbmax,lmax)
      double precision axc(nbmax,lmax),ay2c(nbmax,lmax)
      double precision axcpr(nbmax,lmax),ay2cpr(nbmax,lmax)
      double precision axpr(nbmax,lmax),ay2pr(nbmax,lmax)
      double precision abetac(nbmax,lmax)

      integer l(lmax),ntc(lmax),l0
      integer l0file,lfile(lmax)
      character*500 filename

      common /lvalues1/ l,l0,lmo
      save /lvalues1/


      common /ujlgen1/ abeta,ax,axpr,ay2,ay2pr
      common /ujlgen2/ betamin,dlnbeta
      common /ujlgen3/ abetac,axc,axcpr,ay2c,ay2cpr
      common /ujlgen4/ ntot,ntc,Npoint
      save /ujlgen1/
      save /ujlgen2/
      save /ujlgen3/
      save /ujlgen4/

      open(unit=11,file=filename
     &     ,status='unknown',form='unformatted')
      rewind 11
      read(11)lmofile
      read(11)betamaxfile
      read(11)betamin
      read(11)dlnbeta
      read(11)Npoint

c     TESTING
c      write(123,*)betamaxfile,betamin

      if (lmo.gt.lmofile) then
         write(*,*)'You have entered a lmax'
         write(*,*)'inconsistent with those in the file'
         write(*,*)lmofile
         write(*,*)'You will have to start again'
         stop
      end if

c     Checking if the lvalues.inc file used to build jl file
c     is the same as the one in the code.
        read(11)l0file
        do j=1,l0file
           read(11)lfile(j)
        end do

        do j=1,l0
           if (l(j).ne.lfile(j)) then
              write(*,*)'lvalues.inc file used to build jl file'
              write(*,*)'and the one in the code differ.'
              write(*,*)'You must use the same one'
              write(*,*)j,l(j),lfile(j)
              stop
           end if
        end do

c
c     reading  uj_l information
c     remember to create ujl.dat with ujlgen.f first using
c     correct lmax

c     Number of betas for open models.
      ntot=int(log(betamaxfile/betamin)*dlnbeta)
      do 40 i=1,ntot
         abeta(i)=betamin*exp(dble(i-1)/dlnbeta)
 40   continue

c     Read open models

      do 50 j=1,l0
         do 60 i=1,ntot
            read(11)ax(i,j)
            read(11)ay2(i,j)
 60      continue

c     Closed models

c     The list of betas for closed models will be
c     l+1,l+2,...,l+Npoint,beta(ist+1),beta(ist+2),....,beta(ntot)
c     where ist is given by
         aist=1.0d0+dlnbeta*log(dble(l(j)+Npoint+1)/betamin)
         ist=int(aist)
         ntotc=ntot-ist+Npoint
c     Note that ntotc depends on l. We keep this numbers in an
c     array
         ntc(j)=ntotc

         do 70 i=1,ntotc

c     So that the betas for this l are
            if (i.le.(ntotc-Npoint)) then
               abetac(i,j)=abeta(ntot-i+1)
               abetac(i,j)=dble(int(abetac(i,j)))
            else
               abetac(i,j)=dble(l(j)+1+ntotc-i)
            end if

            read(11)axc(i,j)
            read(11)ay2c(i,j)
 70      continue
 50   continue
      close(11)

c
c     get the interpolation matrices for starting value and deriv.

      do 80 j=1,l0
         call spline(abeta,ax(1,j),ntot,d0lo,d0hi,axpr(1,j))
         call spline(abeta,ay2(1,j),ntot,d0lo,d0hi,ay2pr(1,j))
         call spline(abetac(1,j),axc(1,j),ntc(j),d0lo,
     $        d0hi,axcpr(1,j))
         call spline(abetac(1,j),ay2c(1,j),ntc(j),d0lo,
     $        d0hi,ay2cpr(1,j))
 80   continue

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine powersopen(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the scalar perturbations in the open case.
c Now it is set to a power law in k which is slightly more
c complicated for beta(=ak).

      implicit double precision(a-h,o-z)
      include 'cmbfast.inc'

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn
      common /curvature/ curv,r,k
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr

c     K splitting
      common /cutk/ aksplit, kcutflag
      if (kcutflag.eq.0) then
         win=1
      else
         win=2.0d0*(ak/aksplit)**4
         win=2.0d0*exp(-win)/(1.0d0+exp(-win))
         if (kcutflag.eq.-1) win=1.0d0-win
      end if

c     Normalize so that tilt leaves the same power
c     at k=0.05 Mpc
      anorm=0.05d0
      anorm2=anorm*anorm
      akn=ak/anorm
      aknlog=0.5d0*log((-1.0d0*k/r/r+ak*ak)/anorm2)
      apower=(ak*ak*r*r-k*4.0d0)**2/(-k+ak*r*ak*r)**3
     2                 *r**2*ak*ak*akn*exp((an(in)-2.0d0
     3                  +.5*alphans(in)*aknlog)*aknlog)

c     K splitting
      apower=apower*win

#ifdef OLDNORM
      apower=apower*(4.0d0*3.14159265d0)**2/2.0d0*100.0d0
#endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine powertopen(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the scalar perturbations in the open case.
c Now it is set to a power law in k which is slightly more
c complicated for beta(=ak).

      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      parameter (pi=3.14159265d0)
      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $     itflag,lmaxt
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr
      common /curvature/ curv,r,k

#ifdef OLDNORM
      anorm=0.002d0
#else
      anorm=0.05d0
#endif
      anorm2=anorm*anorm
      akn=ak/anorm

      aux=tanh(pi*ak*r/2.0d0)
      aknlog=0.5d0*log((-3.0d0*k/r/r+ak*ak)/anorm2)
      apower=(ak*ak*r*r-4.0d0*k)/(-k+ak*r*ak*r)
     &*aux/akn*exp((ant(in)+1.0+.5*alphant(in)*aknlog)*aknlog)

#ifdef OLDNORM
      apower=apower*(4.0d0*3.14159265d0)**2*8.0d0*2.0d4
#endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Function Sinhk(x)
      implicit double precision(a-h,o-z)
      common /curvature/ curv,r,k
      if (k.gt.0.0) then
         sinhk=sin(x)
      else
         if (k.lt.0.0) then
c            sinhk=sinh(x)
            y=exp(x)
            sinhk=0.5*(y**2-1.0)/y
         else
            sinhk=x
         endif
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function coshk(x)
      implicit double precision(a-h,o-z)
      common /curvature/ curv,r,k
      if (k.gt.0.0) then
       coshk=cos(x)
      else
       if (k.lt.0.0) then
        coshk=cosh(x)
       else
        coshk=1.0d0
       endif
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initchi(beta,j,ll,chi0,y1,y2,K)
c     For each a given beta and ll=l(j) find the initial
c     chi0 and give y1 and y2 to start the integration
c     of the bessel functions.

      implicit double precision (a-h,o-z)

      include 'cmbfast.inc'

      parameter (nbmax=1000)

      double precision abeta(nbmax),ax(nbmax,lmax)
      double precision ay2(nbmax,lmax)
      double precision axc(nbmax,lmax),ay2c(nbmax,lmax)
      double precision axcpr(nbmax,lmax),ay2cpr(nbmax,lmax)
      double precision axpr(nbmax,lmax),ay2pr(nbmax,lmax)
      double precision abetac(nbmax,lmax)

      integer ntc(lmax)

      common /ujlgen1/ abeta,ax,axpr,ay2,ay2pr
      common /ujlgen2/ betamin,dlnbeta
      common /ujlgen3/ abetac,axc,axcpr,ay2c,ay2cpr
      common /ujlgen4/ ntot,ntc,Npoint

      save /ujlgen1/
      save /ujlgen2/
      save /ujlgen3/
      save /ujlgen4/



      betam1=1.0d0/beta
c     Interpolate initial conditions

c     Open models
      if (K.eq.-1) then


         if ((beta.lt.abeta(ntot)).and.(beta.gt.abeta(1))) then

c     Beta is inside array
            ibl=int(dlnbeta*log(beta/betamin))
            ibh=ibl+1
            ho=abeta(ibh)-abeta(ibl)
            a0=(abeta(ibh)-beta)/ho
            b0=(beta-abeta(ibl))/ho
         else
c     Extrapolate with a constant
            if (beta.gt.abeta(ntot)) then
               ibl=ntot
               ibh=ntot
               a0=1.0d0
               b0=0.0d0
               ho=0.0d0
            else
c     Extrapolate linearly
               ibl=1
               ibh=ibl+1
               ho=abeta(ibh)-abeta(ibl)
               a0=(abeta(ibh)-beta)/ho
               b0=(beta-abeta(ibl))/ho
               ho=0.0d0
            end if
         end if
         aux2=a0*ax(ibl,j)+b0*ax(ibh,j)+((a0**3-a0)*axpr(ibl,j)
     &        +(b0**3-b0)*axpr(ibh,j))*ho*ho/6.d0
         aux1=exp(aux2)*betam1
         chi0=log(aux1+sqrt(aux1**2+1.0d0))
         y1=1.0d-6*aux1/dble(ll)
         y2=a0*ay2(ibl,j)+b0*ay2(ibh,j)+((a0**3-a0)*ay2pr(ibl,j)
     &        +(b0**3-b0)*ay2pr(ibh,j))*ho*ho/6.d0

      else
         if (beta.lt.abetac(ntc(j)-Npoint+1,j)) then
c     Beta is inside first part of array
            ibh=ll+1+ntc(j)-int(beta)
            ibl=ibh-1
            ho=abetac(ibh,j)-abetac(ibl,j)
            a0=(abetac(ibh,j)-beta)/ho
            b0=(beta-abetac(ibl,j))/ho
         else
            if (beta.lt.abetac(1,j)) then
c     Beta is inside second part of array
               ib1=int(dlnbeta*log(beta/betamin))
               ibh=ntot-ib1
               ibl=ibh-1
               if (beta.gt.abetac(ibl,j)) then
                  ibl=ibl-1
                  ibh=ibh-1
               end if
               ho=abetac(ibh,j)-abetac(ibl,j)
               a0=(abetac(ibh,j)-beta)/ho
               b0=(beta-abetac(ibl,j))/ho
            else
               ibh=1
               ibl=1
               a0=1.0d0
               b0=0.0d0
               ho=0.0d0
            end if
         end if

         aux2=a0*axc(ibl,j)+b0*axc(ibh,j)+((a0**3-a0)*axcpr(ibl,j)
     &        +(b0**3-b0)*axcpr(ibh,j))*ho*ho/6.d0
         aux1=exp(aux2)*betam1
         chi0=asin(aux1)
         y1=1.0d-6*aux1/dble(ll)
         y2=a0*ay2c(ibl,j)+b0*ay2c(ibh,j)+((a0**3-a0)*ay2cpr(ibl,j)
     &        +(b0**3-b0)*ay2cpr(ibh,j))*ho*ho/6.d0

      end if



      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine evalujl(chi,y1,y2,delchi,beta2,ap1)

c     Given y1 and y2, at chi, the function and its
c     derivative, move forward one step to chi+delchi.
c     At output chi, y1, y2 will be updated.
c     K=-1 for open and +1 for closed. ap1=l*(l+1)
c     and beta2 is beta**2

      implicit double precision (a-h,o-z)


      dydchi1=y2
      dydchi2=(ap1/sinhk(chi)**2 - beta2)*y1
      hh=delchi*0.5d0
      h6=delchi/6.0d0
      xh=chi+hh
      yt1=y1+hh*dydchi1
      yt2=y2+hh*dydchi2
      dyt1=yt2
      dyt2=(ap1/sinhk(xh)**2 - beta2)*yt1
      yt1=y1+hh*dyt1
      yt2=y2+hh*dyt2
      dym1=yt2
      dym2=(ap1/sinhk(xh)**2 - beta2)*yt1
      yt1=y1+delchi*dym1
      dym1=dyt1+dym1
      yt2=y2+delchi*dym2
      dym2=dyt2+dym2
      dyt1=yt2
      dyt2=(ap1/sinhk(chi+delchi)**2 - beta2)*yt1
      y1=y1+h6*(dydchi1+dyt1+2.0d0*dym1)
      y2=y2+h6*(dydchi2+dyt2+2.0d0*dym2)

      chi=chi+delchi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine symsource(s,iflip,nstep)
      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

c      parameter (nstep0=2400)
      double precision s(nstep0,2)
      integer iflip

      do i=1,iflip-1
         stemp1=s(2*iflip-i,1)
         stemp2=s(i,1)
         s(2*iflip-i,1)=stemp1+stemp2
         s(2*iflip-i,2)=stemp1-stemp2
         s(i,1)=s(2*iflip-i,1)
         s(i,2)=-s(2*iflip-i,2)
      end do
      s(iflip,1)=2.0d0*s(iflip,1)
      s(iflip,2)=0.0d0
      do i=2*iflip,nstep
         s(i,2)=s(i,1)
      end do
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine symtau(r,atau0,dtaureg,nreg,nr,irec,nstep0)

c     Symmetrizes the time array around chi=pi/2
c     It can only take 10 subdivisions.

      implicit double precision(a-h,o-z)

      double precision atau0(*),dtaureg(*)
      integer nreg(*),nr,nstep0
      parameter(pi=3.14159265359d0)
      parameter(pio2=pi/2.0d0)
c     Internal variables

      double precision dtauf(20),dtaub(20),chi(20),tau(20)
      integer nd,ipoint(20)

      common /tnowkb/twkb1,twkb2,twkb3,twkb4

c     TESTING
c      write(70,*)'r=',r

      nd=nr+1
      tau0=atau0(nreg(nd))
      taustart=atau0(1)

      if (r.le.0.0d0) return
      chi0=(tau0-atau0(nreg(1)))/r
      if (chi0.le.pio2) return

c     We do not want to switch to wkb inside recombination. Remember
c     times to check later.
      chiwkb1=(tau0-twkb1)/r
      chiwkb2=(tau0-twkb2)/r
      if (chiwkb1.gt.pio2) then
         chiwkb1=pi-chiwkb1
         wkbtemp1=tau0-chiwkb1*r
      end if
      if (chiwkb2.gt.pio2) then
         chiwkb2=pi-chiwkb2
         wkbtemp2=tau0-chiwkb2*r
      end if
      twkb3=min(wkbtemp1,wkbtemp2)
      twkb4=max(wkbtemp1,wkbtemp2)


c     TESTING
c         write(70,*)'chi,tau,dtauf,dtaub'

      do i=1,nd
         ipoint(i)=i
         chi(i)=(tau0-atau0(nreg(i)))/r
         tau(i)=atau0(nreg(i))
         j=min(i,nr)
         dtauf(i)=dtaureg(j)
         j=mod(i-2+nr,nr)+1
         dtaub(i)=dtaureg(j)
c
         if (chi(i).ge.pi) then
            write(*,*)'Parameters not allowed in symtau'
            write(*,*)'Inside region of death'
            write(*,*)'(tau0,tau,r,n)='
            write(*,*)tau0,atau0(nreg(i)),r,nreg(i)
            stop
         end if
c     TESTING
c         write(70,'(4E15.5)')chi(i),tau(i),dtauf(i),dtaub(i)

      end do

      id=0
      do i=1,nd
         if (chi(i).gt.pio2) then
            id=i
            chi(i)=pi-chi(i)
            tau(i)=tau0-chi(i)*r
            dtau=dtauf(i)
            dtauf(i)=dtaub(i)
            dtaub(i)=dtau
         end if
      end do

c     TESTING
c      write(70,*)'id=',id
c      write(70,*)'A,chi,tau,dtauf,dtaub'
c      do i=1,nd
c         write(70,'(4E15.5,1I6)')chi(i),tau(i),dtauf(i)
c     $        ,dtaub(i),ipoint(i)
c      end do

      do i=1,id
         if (tau(i).lt.tau(id+1)) then
            dtauf(i)=min(dtauf(i),dtaub(id+1))
            dtaub(i)=min(dtaub(i),dtaub(id+1))
         end if
         do j=id+1,nd-1
            if ((tau(i).ge.tau(j)).and.(tau(i).lt.tau(j+1))) then
               dtauf(i)=min(dtauf(i),dtaub(j+1))
               dtaub(i)=min(dtaub(i),dtauf(j))
            end if
         end do
      end do

c     TESTING
c      write(70,*)'B,chi,tau,dtauf,dtaub'
c      do i=1,nd
c         write(70,'(4E15.5,1I6)')chi(i),tau(i),dtauf(i)
c     $        ,dtaub(i),ipoint(i)
c      end do


      nr=nr+1
      nd=nd+1
      chi(nd)=pio2
      tau(nd)=tau0-pio2*r
      dtauf(nd)=dtaub(nd-1)
      dtaub(nd)=dtaub(nd-1)
      ipoint(nd)=nd

c     TESTING
c      write(70,*)'chi,tau,dtauf,dtaub'
c      do i=1,nd
c         write(50,'(4E15.5,1I6)')chi(i),tau(i),dtauf(i)
c     $        ,dtaub(i),ipoint(i)
c      end do

c     From the array of tau, construct region arrays

      call tau2reg(tau,dtauf,dtaub,dtaureg,nreg
     $     ,ipoint,nr)
      nd=nr+1

c     TESTING
c      write(70,*)'chi,tau,dtauf,dtaub'
c      do i=1,nd
c         write(70,'(4E15.5,1I6)')chi(i),tau(i),dtauf(i)
c     $        ,dtaub(i),ipoint(i)
c      end do

c     TESTING
c      write(70,*)'nreg'
c      do i=1,nd
c         write(70,*)i,nreg(i),ipoint(i)
c      end do
c     TESTING
c      write(70,*)'dtaureg'
c      do i=1,nr
c         write(70,*)i,dtaureg(i)
c      end do
c
c      call n2tau(atau0,tau(nd),dtaureg,nreg,nr,nstep0)
c      nd=nr+1
c
c     TESTING
c      write(70,*)'taus'
c      do i=1,nd
c         write(70,*)nreg(i),atau0(nreg(i))
c      end do

c     Constructing the whole tau.

c     Find where the time corresponding to
c     the start of recombination is.
      irec=0
      ifl=1
      do while (irec.eq.0)
         do i=1,nd
            if (ipoint(i).eq.ifl) then
               irec=i
               goto 20
            end if
         end do
 20      continue
         ifl=ifl+1
      end do

c     TESTING
c      write(70,*)'irec',irec

      nr=nr+irec-1
      nd=nd+irec-1

      do j=nd,irec,-1
         nreg(j)=nreg(j-irec+1)+nreg(irec)-1
      end do
      do j=irec-1,1,-1
         nreg(j)=nreg(2*irec-1)-nreg(2*irec-j)+1
      end do

c     TESTING
c      write(70,*)'nreg(irec)',nreg(irec)
c
c      write(70,*)'nreg sym'
c      do i=1,nd
c         write(70,*)i,nreg(i)
c      end do

c     TESTING
c      write(70,*)'j,int(abs(real(j-irec)+0.5))+irec'
      do j=nr,irec,-1
         i=int(abs(real(j-irec)+0.5))+1
         dtaureg(j)=dtaureg(i)
c     TESTING
c         write(70,*)j,i
      end do
      do j=irec-1,1,-1
         i=int(abs(real(j-irec)+0.5))+irec
         dtaureg(j)=dtaureg(i)
c     TESTING
c         write(70,*)j,i
      end do

c     TESTING
c      write(70,*)'dtaureg sym'
c      do i=1,nr
c         write(70,*)i,dtaureg(i)
c      end do

c     TESTING
c      write(70,*)'taustart,nr,nstep0'
c      write(70,*)taustart,nr,nstep0

      call n2tau(atau0,taustart,dtaureg,nreg,nr,nstep0)
c     TESTING
c      write(70,*)'taus'
c      do i=1,nd
c         write(70,*)nreg(i),atau0(nreg(i))
c      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tau2reg(tau,dtauf,dtaub,dtaureg,nreg
     $     ,ipoint,nr)

      implicit double precision(a-h,o-z)
      double precision tau(*),dtauf(*),dtaub(*),dtaureg(*)
      integer nreg(*),nr
      integer ipoint(*)

      nd=nr+1
      do i=1,nd
         ipoint(i)=i
      end do

c     TESTING
c      write(70,*)'ipoint,tau(ipoint(i)),dtauf,dtaub'
c      do i=1,nd
c         write(70,'(1I6,3E15.5)')ipoint(i),tau(ipoint(i))
c     $        ,dtauf(ipoint(i)),dtaub(ipoint(i))
c      end do

c     Ordering array with respect to tau,  dum way
      do j=2,nd
         tautemp=tau(ipoint(j))
         iptemp=ipoint(j)
         do i=j-1,1,-1
            if (tau(ipoint(i)).le.tautemp) goto 10
            ipoint(i+1)=ipoint(i)
         end do
         i=0
 10      ipoint(i+1)=iptemp
      end do

c     TESTING
c      write(70,*)'ipoint,tau(ipoint(i)),dtauf,dtaub'
c      do i=1,nd
c         write(70,'(1I6,3E15.5)')ipoint(i),tau(ipoint(i))
c     $        ,dtauf(ipoint(i)),dtaub(ipoint(i))
c      end do

c     Eliminate unnecessary regions
      i=2
      do while (i.le.nd-1)
         deltatau=tau(ipoint(i))-tau(ipoint(i-1))
         dtau=min(dtauf(ipoint(i-1)),dtaub(ipoint(i)))
         ntemp=int(deltatau/dtau+0.5d0)
         dtau2=min(dtauf(ipoint(i)),dtaub(ipoint(i+1)))
         reldtau=abs(dtau/dtau2-1)
         if (((ntemp.lt.5).or.(reldtau.lt.0.1)).and.
     $        (ipoint(i).ne.1)) then
            dtauf(ipoint(i-1))=min(dtauf(ipoint(i))
     $           ,dtauf(ipoint(i-1)))
            dtaub(ipoint(i-1))=min(dtaub(ipoint(i))
     $           ,dtaub(ipoint(i-1)))
            do j=i,nd-1
               ipoint(j)=ipoint(j+1)
            end do
            nd=nd-1
            i=i-1
         end if
         i=i+1
      end do
      nr=nd-1

c     TESTING
c      write(70,*)'ipoint,tau(ipoint(i)),dtauf,dtaub'
c      do i=1,nd
c         write(70,'(1I6,3E15.5)')ipoint(i),tau(ipoint(i))
c     $        ,dtauf(ipoint(i)),dtaub(ipoint(i))
c      end do

c     Constucting dtaureg, nreg
      nreg(1)=1
      do i=1,nr

c     Initial and final times for region
         tauinitial=tau(ipoint(i))
         taufinal=tau(ipoint(i+1))

c     Timestep for region
         dtau=min(dtauf(ipoint(i)),dtaub(ipoint(i+1)))
         if (dtau.le.0) then
            write(*,*)'Error in symtau, dtau<0'
            stop
         end if
         an=(taufinal-tauinitial)/dtau+0.5d0
         if (an.lt.0.0d0) then
            write(*,*)'Error in symtau, an<0'
            stop
         end if
         n=int(an)
         n=max(n,1)
         dtau=(taufinal-tauinitial)/dble(n)

c     Constructing region and step arrays
         dtaureg(i)=dtau
         nreg(i+1)=nreg(i)+n
      end do

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tau2n(tau,n,atau0,dtaureg,nreg,nr)

c     Given a tau finds the position in the atau0 array

      implicit double precision(a-h,o-z)

      double precision atau0(*),dtaureg(*)
      integer nreg(*),nr

      do i=1,nr
         if ((tau.lt.atau0(nreg(i+1))).and.
     $        (tau.ge.atau0(nreg(i)))) then
            an=(tau-atau0(nreg(i)))/dtaureg(i)
            n=nreg(i)+int(an)
            goto 10
         end if
      end do
      if (tau.eq.atau0(nreg(nr+1))) then
         n=nreg(nr+1)
      else
         write(*,*)'tau out of range in tau2n'
         write(*,*)tau,atau0(nreg(nr+1))
         stop
      end if
 10   continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine n2tau(atau0,taustart,dtaureg,nreg,nr,nstep0)

c     Creates a atau0 array out of the region information.

      implicit double precision(a-h,o-z)

      double precision atau0(*),dtaureg(*)
      integer nreg(*),nr,nstep0

      nd=nr+1
      nstep=nreg(nd)
      if (nstep.gt.nstep0) then
         write(*,*)
     2        'Sorry, the arrays were dimensioned for a max. of',
     3        nstep0, 'timesteps.',
     4        'The model you requested needs',
     5        nstep, 'Please make the arrays bigger by making ',
     6        'nstep0 bigger where it appears'
         stop
      end if

      atau0(1)=taustart
      do i=1,nr
         do j=nreg(i)+1,nreg(i+1)
            atau0(j)=atau0(j-1)+dtaureg(i)
         end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
