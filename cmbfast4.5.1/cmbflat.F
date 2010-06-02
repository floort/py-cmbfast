C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (mzaldarriaga@cfa.harvard.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cmbflat(clts,cltt,cles,
     &     clet,clbt,clcs,clct,clkk,cltk)

c     This program integrates the linearized equations of general relativity,
c     the Boltzmann equations and the fluid equations for scalar perturbations
c     of a spatially flat Friedmann-Robertson-Walker universe in synchronous
c     gauge.  The time variable is conformal time dtau=dt/a(t) and the spatial
c     dependence is Fourier transformed.  ak is the comoving wavevector;
c     comoving distances are x=r/a(t), with a(t)=1 today.  The units of both
c     length and time are Mpc.
c
c
c     All parameters are passed in common statements.
      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'
#ifdef MPIBUILD
      include 'mpif.h'   !MPI include statements
#endif



      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr
      common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu
     &     ,nqmax,iq0,iq1,iq2
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     $     grhonr
      common /timesteps/ atau0,dtau1,dtau2
      common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
      common /reionization2/ j2ri1,nri,nri0

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn
      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $       itflag,lmaxt

c     Dependent integration variables: a,phi,delta_c,theta_c,delta_b,theta_b,
c     2*(lmaxg+1) photon moments (total intensity plus polarization),
c     (lmaxnr+1) massless neutrino moments, (lmaxnu+1) massive neutrino
c     moments, nqmax momenta for each l.
c
c     For the tensor modes me will calculate a, ht, htprime, and
c     2*(lmaxt+1) photon moments (total intensity plus polarization).
c
c     The metric perturbation variables are computed algebraically.
      parameter (lmax0=12,lmaxnr0=25,lmaxnu0=25,nqmax0=15,lmaxt0=10)
      parameter (lm1=8,lm2=7,lm3=4,nq1=15)
c     lmx0 MUST be grater than all of the previous lmax0.
      parameter (lmx0=30)

      common /qparm/ wdyn,ndyn
      double precision dynrho,wdyn_func
      external dynrho,wdyn_func,dyn_phi,dyn_nrg

      parameter (nvar0=9+2*(lmax0+1)+(lmaxnr0+1)+nqmax0*(lmaxnu0+1))

      parameter (nvar0t=2*(lmaxt0+1)+2+1+lmaxnr0+1)
      parameter (ketamax0=3*l0max+126)
      parameter (nkmax=2*l0max)


c     Transfer functions

      double precision tautf(ntfmax),ztf(ntfmax)
      double precision tautf0(ntau+ntfmax)
      integer indtau(ntau+ntfmax)
      common /lens/ lensflag
      common /transfer/akmaxt,ztf,nlnkt,ict,ntf
      double precision power(nk0*6,ntfmax),phi,phi1,phihk
      integer istore(ntfmax)
      common/transfer_function/power,istore
      double precision atrg(ntau,4),akl(nstep0),phik(nstep0,ntau)
c      common /par/atrg,taur,tau0,akl,phik,nkt

c
c     Scalar perturbations
      double precision y(nvar0),yprime(nvar0)
      double precision d(nk0,nstep0),dpr(nk0,nstep0)
      double precision dp(nk0,nstep0),dppr(nk0,nstep0)
      double precision dkk(nk0,nstep0),dkpr(nk0,nstep0)
      double precision dl(lmax),dl2(lmax),dl3(lmax)
      double precision dpl(lmax),dpl2(lmax),dpl3(lmax)
      double precision dkl(lmax),dkl2(lmax),dkl3(lmax)
      double precision s2(nstep0),sp2(nstep0),sk2(nstep0)
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
      double precision dtl(lmax),dtl2(lmax),dtl3(lmax)
      double precision dtel(lmax),dtel2(lmax),dtel3(lmax)
      double precision dtbl(lmax),dtbl2(lmax),dtbl3(lmax)
      double precision st2(nstep0),ste2(nstep0)
      double precision stb2(nstep0)
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

c
      double precision atau0(nstep0),dtau1(nstep0),dtau2(nstep0)
      double precision ak0(nk0)
      double precision ajl(ketamax0,lmax),ajlpr(ketamax0,lmax)
      double precision xx(ketamax0), dxx(ketamax0)
      integer mxx(nstep0)
      double precision ak1(nkmax),dak1(nkmax)
      parameter (d0hi=1.0d40,d0lo=1.0d40)
      integer l(lmax)
      double precision xl(lmax),temp,growthfactor,akgf

      common /lvalues1/ l,l0,lmo
      common /lvalues2/ akmax0
      save /lvalues1/
      save /lvalues2/

      common /jlgen/ ajl,ajlpr,xx,dxx,kmax
      save /jlgen/

c     TESTING. This arrays can be used output the contributions
c     of different ks to the cls.
c     double precision ctt(lmax0),cpp(lmax0)

c
      common /epsilon/ epsw
c
c
c     const=7*pi**4/120.
      parameter (const=5.68219698d0,zeta3=1.20205690d0)
c
      parameter (fourpi=4.0d0*3.14159265d0,xlimmin=10.0d0)
c     dverk  parameters.
      parameter (tol=1.0d-8,nw=nvar0,nwt=nvar0t)
      dimension c(24),w(nw,9),wt(nwt,9)
c
c
      dimension denl(lmx0),dlfdlq(nqmax0)
      common /store/ denl,dlfdlq

      common /maxreion/ taurmax,armax,jrmax
c
      double precision dtauda, dsoundda
      external fderivs,fderivst,dtauda,dsoundda

cxxxc     TIMING
cxxxc     Timing variables for testing purposes
cxxxc     Should be commented in released versions
cxxx      real etime,actual,timeprev,tarray(2)

c     TO SAVE MEMORY
      common /memory1/ d,dpr,dp,dppr,dkk,dkpr
      common /memory2/dt,dtpr,dte,dtepr,dtb,dtbpr

c     K splitting parameters
      common /cutk/ aksplit, kcutflag

cxxxc     TIMING
cxxxc     Enterning CMBFLAT
cxxx      actual=etime(tarray)

#ifdef MPIBUILD
        CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )
        CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
#else
        myid=0
        numprocs=1
#endif
c      if (lensflag.ne.0) write(17,*)ntau,h0/100.0
      omegam=omegab+omegac+omegan
c     zst is the value of z where the program stops the calculation.
      if ((abs(omegab+omegac-1.0d0).lt.1.0d-3).and.(zri.eq.0.0).
     &     and.(optdlss.eq.0).and.(itflag.eq.0)
     $     .and.(h0.gt.40.0)) then
         zst=10.0d0
      else
         zst=0.0d0
      endif
      if (ict.ne.0) zst=min(zst,ztf(ntf))
c
      if (itflag.ne.2) then
         lmaxg=lmax0
      end if
#ifdef HP
c     At l=2000 increasing lmaxnr even further from 25 to 50 makes
c     a small difference, less than 0.1%
c     use same lmaxnr for scalars and tensors
            lmaxnr=lmaxnr0
#else
            lmaxnr=lm2
#endif
      nstep=nstep0
      if (itflag.ne.0) then
         lmaxt=lmaxt0
      else
         lmaxt=0
      end if
c

      do 133 j=1,lmx0
         denl(j)=1.0d0/dble((2*j+1))
 133  end do
c
c
c     Initialize neutrino mass and storage.
      if (annunr.eq.0.or.omegan.eq.0.0) then
         amnu=0.0d0
         nqmax=0
         lmaxnu=0
      else
c     amnu=m_nu*c**2/(k_B*T_nu0).
c     The 1.0e-18 is to prevent certain compilers (e.g. Sun's) from
c     thinking that a divide by zero will occur.
      amnu=const/(1.5d0*zeta3)*grhom/grhonr*omegan
     &        /dble((annunr+1.0d-18))

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
         nvart=(2*lmaxt+5)+lmaxnr+1
      else
         nvart=0
      end if

c     Initialize massive neutrinos.
      if ((annunr.ne.0).or.(omegan.ne.0.0)) then
         arel=1.d-3/amnu
         call initnu1(amnu)
      end if

c     Time today
      tau0=rombint(dtauda,1.0d-8,1.0d0,tol)
      epsw=100.0d0/tau0
c Calculating when reinization starts.
      if (optdlss.gt.0.0d0) then
         call reiopar(optdlss,zri,zristp,rif)
      end if

c     TESTING
c      write(*,*)optdlss,zri,zristp,rif

c     Distance to LSS
        dlss=tau0

c
c
c     Maximum and minimum k-values. dlnk is the logarithmic spacing
c     used for low k.
      akmax=akmax0/tau0
      akmin=0.15d0/tau0

c     TESTING
c      akmin=1.0d-4
c      akmax=0.22
c      write(*,*)akmax,tau0

      if (itflag.eq.0) then
         dlnk=0.1d0
      else
         dlnk=0.05d0
      endif
c
c
c     Timesteps during recombination (tentative, the actual
c     timestep is the minimum between this value and taurst/40,
c     where taurst is the time when recombination starts.
c     akmax is the maximum ak being calculated, 4.0/akmax
c     seems a very poor sampling, but akmax is bigger than
c     lmax/dlss, this is why this condition works.

      if (akmax.ne.0) then
         dtaurec=4.0d0/akmax
      else
         dtaurec=0.0d0
      end if
c     Stoping time
      taumax=rombint(dtauda,1.0d-8,1.0d0/(zst+1.0d0),tol)

c     TESTING
c      write(*,*)'taumax=',taumax,zst

c     Timesteps after recombination, exponentially separated.

      dlntau0=0.501d-2

      if (itflag.ne.0) dlntau0=0.025d-1
c     Time at which reionization takes place
      if (zri.ne.0.0) then
         taurist=rombint(dtauda,1.0d-8,1.0d0/(1+zri),tol)
         tauristp=rombint(dtauda,1.0d-8,1.0d0/(1+zristp),tol)
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
      call finithermo(taumin,taumax,tau0,taurend,dlntau0,n1,nstep)

c     TESTING
c      write(*,*)'tau0,atau0(nstep)'
c      write(*,*)tau0,atau0(nstep)

      dlsso=tau0-taurmax
c     Sound Horizon at LSS
      shor=rombint(dsoundda,1.0d-8,armax,tol)

c     K splitting
c     Compute cut wavevector
      aksplit=1.5d0/shor*aksplit
#ifdef MAX

      open(unit=28,file='dlss.dat',
     &     status='unknown',form='formatted')
      rewind 28
      write(28,'(6E13.5)')dlsso,taurmax,tau0,shor,armax,aksplit
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


c
c     Calculating the times for the outputs of the transfer functions.
c
      if (ict.ne.0) then
            a0=1.0d-8
         do 150 itf=1,ntf
            atf=1.0d0/(ztf(itf)+1.0d0)
            tautf(itf)=rombint(dtauda,a0,atf,tol)
            tautf(itf)=min(tautf(itf),atau0(nstep))
            indtau(itf)=itf
 150     enddo
         if (lensflag.ne.0) then
            do itf=1,ntf
             tautf0(itf)=tautf(itf)
            enddo
            do i=1,ntau+ntf
             indtau(i)=0
            enddo
            ar=1.0d0/1090.0d0
            taur=rombint(dtauda,a0,ar,tol)
            itf=1
            do i=1,ntau
              taunew=dble(i-1)*(tau0-taur)/dble(ntau-1)+taur
112         continue
              if (taunew.le.tautf0(itf).or.itf.gt.ntf) then
                tautf(i+itf-1)=taunew
              else
                tautf(i+itf-1)=tautf0(itf)
                indtau(i+itf-1)=itf
                itf=itf+1
                goto 112
              endif
            enddo
            if (itf.le.ntf) then 
             do i=itf,ntf
              indtau(ntau+i)=i
              tautf(ntau+i)=tautf0(i)
             enddo
            endif
            ntf=ntf+ntau
         endif
      endif
      tautf(ntf-1)=tautf(ntf-1)-tol       
c
c     Integration will only be carried out after z=10 for low
c     k, k < k10. If there is reionization the boundary will not be
c     z=10 but tauristp
      if (zst.ne.10.0d0) then
         t10=rombint(dtauda,1.0d-8,1.0d0/11.0d0,tol)
c         ak10=250.0d0/tau0
c     TESTING
         ak10=500.0d0/tau0

         n10=n1+int(log(t10/taurend)/dlntau0)
         if (zri.ne.0.0) then
            nriend=j2ri1+nri+nri0
            t10=max(t10,atau0(nriend))
            n10=nriend+int(log(t10/atau0(nriend))/dlntau0)
         end if
c     TESTING
c         write(*,*)'Increasing ak10'
c         ak10=ak10*2.0
c         Do not stop
c         t10=tau0
c         ak10=akmax
c         n10=nstep+1
      else
         t10=tau0
         ak10=akmax
         n10=nstep+1
      end if


cxxxc     TIMING
cxxx      timeprev=actual
cxxx      actual=etime(tarray)
cxxx      write(*,*)actual-timeprev,' END INITIALIZE, NOW SOURCES'
cxxx      write(50,*)actual-timeprev,' END INITIALIZE, NOW SOURCES'
cxxx
c
c     Calculation of the CMB sources.
c
      nk=0
      if (ict.ne.1) then
c     set k values for which the sources for the anisotropy and
c     polarization will be calculated. For low values of k we
c     use a logarithmic spacing.
c
         if (zri.ne.0.0) then
#ifdef HP
            dlnk0=1.0d0*dlnk
#else
            dlnk0=2.0d0*dlnk
#endif
         else
            dlnk0=5.0d0*dlnk
         end if
#ifdef HP
         dkn1=0.4d0/taurmax
         dkn2=0.8d0/taurmax


#else
         dkn1=0.8d0/taurmax
         dkn2=1.5d0/taurmax
#endif

         nk1=int(log(dkn1/akmin/dlnk0)/dlnk0)+1

cxxx         if (akmax.gt.(1500.0d0/tau0)) then
cxxx            nk2=int((1500.0d0/tau0-akmin*exp((nk1-1)*dlnk0))
cxxx     &           /dkn1)+nk1+1
cxxx            nk=int((akmax-1500.0d0/tau0)/dkn2)+nk2+1
cxxx         else
cxxx            nk=int((akmax-akmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1
cxxx            nk2=nk+1
cxxx         end if

cxxx         akchange=1500.0d0/tau0
         akchange=5.0d0*3.14159265d0/shor

         if (akmax.gt.akchange) then
            nk2=int((akchange-akmin*exp((nk1-1)*dlnk0))
     &           /dkn1)+nk1+1
            nk=int((akmax-akchange)/dkn2)+nk2+1
         else
            nk=int((akmax-akmin*exp((nk1-1)*dlnk0))/dkn1)+nk1+1
            nk2=nk+1
         end if

         if (nk.gt.nk0) then
            write(*,*)
     2           'Sorry, the arrays were dimensioned for a maximum of'
            write(*,*) nk0, 'k modes.'
            write(*,*)'The model you requested needs',nk
            write(*,*)'Please make the arrays bigger by making '
            write(*,*)'nk0 bigger where it appears'
            stop
         end if

cxxxc     TIMING
cxxx        write(50,*)'ak sources',nk1,nk2,nk
cxxx        write(*,*)'ak sources',nk1,nk2,nk

         do 25 ik=1,nk
            if (ik.le.nk1) then
               ak0(ik)=akmin*exp((ik-1)*dlnk0)
            else
               if (ik.gt.nk2) then
                  ak0(ik)=ak0(nk2)+dble(ik-nk2)*dkn2
               else
                  ak0(ik)=ak0(nk1)+dble(ik-nk1)*dkn1
               end if
            end if
 25      continue

c     TESTING
c     To output average of the sources
c         open(unit=25,file='average_1.dat',status='unknown'
c     $        ,form='formatted')
c         rewind 25

         kkmin=1
         kkmax=nk

c     K split
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
cxxx 1234    continue
cxxx         open(unit=55,file='traces.dat',status='unknown'
cxxx     $        ,form='formatted')
cxxx         rewind 55
cxxx         write(*,*)'Enter ak,tol1,epsilon1,epsilon2'
cxxx         read(*,*)ak,tol1,epsilon1,epsilon2
cxxx         ak2=ak*ak
cxxx         taustart=0.001d0/ak
cxxx         taustart=min(taustart,0.1d0)
cxxx         taustart0=taustart
cxxx         call finitial(y,taustart)
cxxx         tau=taustart
cxxx         ind=1
cxxxc     Begin timestep loop.
cxxx         ntimes=1000
cxxx         tol1=dmin1(tol,1.0e-17*(ak/1.0e-5)**2*tol1)
cxxx         write(*,*)tol1
cxxx         ind=2
cxxx         do j=1,24
cxxx            c(j)=0.0d0
cxxx         end do
cxxx         c(3)=1.0d-8
cxxx         do j=2,ntimes
cxxx            tauend=taustart*
cxxx     $           exp(dble(j-1)*log(tau0/taustart)/dble(ntimes))
cxxxc            tol1=tol
cxxx            call dverk(nvar,fderivs,tau,y,
cxxx     $           tauend,tol1,ind,c,nw,w)
cxxx            write(55,'(20E15.5)')tau,y(4),y(5),y(6),y(7),y(8),y(9)
cxxx     $           ,y(10)
cxxx          end do
cxxx          close(55)
cxxx          goto 1234
cxxxc     END TESTING

c     Loop over wavenumbers.
         do 30 ik=kkmin+myid,kkmax,numprocs  !MPI modification

            ak=ak0(ik)

c     TESTING. Output ks and betas where the source is calculated.
c            write(68,*)ik,ak

c
c     Begin when wave is far outside horizon.
c     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
c     of relativistic neutrinos.
            taustart=0.001d0/ak
c     Make sure to start early in the radiation era. We are not including
c     neutrinos as sources of the tensor modes.
            taustart=min(taustart,0.1d0)
            taustart0=taustart
c     Start when massive neutrinos are strongly relativistic.
            if (amnu.ne.0.0d0) then
               arel=1.d-3/amnu
               taustart=min(taustart0,arel/adotrad)
            end if
c
c     TESTING this variables can be used to compute
c     the sources averaged over time.
c            averaget=0.0d0
c            averagep=0.0d0

            ak2=ak*ak
c

            if (itflag.ne.2) then
               call finitial(y,taustart)

       if(ndyn.eq.1.or.ndyn.eq.2) then
        y(nvar-1) = y(nvar0-1)
        y(nvar)   = y(nvar0)
       endif

               tau=taustart
               ind=1
c     Begin timestep loop.
c
c     d contains the sources for the anisotropy and dp
c     for the polarization. t means tensor.

               if (ict.ne.0) itf=1
               if (lensflag.ne.0) itau=1

c     TESTING
c               tol1=dmin1(tol,1.0e-17*(ak/1.0e-5)**2)
c               write(*,*)ak,tol1

c     The instability reported by David Spergel ocurring at low k was
c     cured by significantly increasing the tolerance. This could lead to
c     complaints by dverk so hmin had to be lowered.

c     TESTING
cxxx               ind=2
cxxx               tol1=tol
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

               do 10 j=2,nstep
                  tauend=atau0(j)
                  if (ak.gt.ak10.and.tauend.gt.t10.and.
     2                 (ict.eq.0.or.tau.gt.tautf(ntf))) then
                     d(ik,j)=0.0d0
                     dp(ik,j)=0.0d0
                  else
                     call dverk(nvar,fderivs,tau,y,tauend,tol1
     &                    ,ind,c,nw,w)
c     TESTING. Output some components of the y vector
c     for a particular k
c
c                     if (ik.eq.2) then
c                        write(69,'(3E15.5)')tau,y(7),y(9)
c                     end if

                     call fderivs(nvar,tau,y,yprime)
                     call foutput(nvar,y,yprime,j,tau0,tau
     &                    ,d(ik,j),dp(ik,j),dkk(ik,j),phi)

cxxxc     TESTING. Computing average sources.
cxxx                     ip1=min(j+1,nstep)
cxxx                     im1=max(1,j-1)
cxxx                     tdtau=(atau0(ip1)-atau0(im1))/dble(ip1-im1)
cxxx                     averaget=averaget+d(ik,j)*tdtau
cxxx                     averagep=averagep+dp(ik,j)*tdtau
c     TESTING. Output of sources at maximum of visibility function.
c     if (j.eq.jrmax) then
c     write(111,*)ak0(ik),d(ik,j),dp(ik,j)
c     end if

c     Calculation of transfer functions.
 101                 continue
                     if (ict.ne.0.and.itf.le.ntf) then
                        if (j.lt.nstep.and.tauend.lt.tautf(itf).
     2                       and.atau0(j+1).gt.tautf(itf)) then
                           call dverk(nvar,fderivs,tau,y,tautf(itf)
     2                          ,tol1,ind,c,nw,w)
                        endif
c     output transfer functions for this k-value.
                   if (abs(tau-tautf(itf)).lt.1.0d-5) then
                     if(indtau(itf).ne.0) then
                      call outtransf(nvar,y,0.0d0,ik,indtau(itf))
                     else
                      if (lensflag.ne.0) then
                       call fderivs(nvar,tau,y,yprime)
                       call foutput(nvar,y,yprime,j,tau0,tau
     &                 ,temp,temp,temp,phi)
                       if (itau.eq.1) akl(ik)=ak
                       if (ik.eq.1) then 
                        atrg(itau,1)=tau
                        atrg(itau,2)=y(1)
                        atrg(itau,3)=tau
                       endif
                       phik(ik,itau)=phi
c                       write(17,*)y(1),tau,phi
c              if (itau.eq.ntau) then 
c                   do in=1,nn
c                     call powersflat(ak,in,apowers)
c                     write(21,'(5E14.5)')ak/(h0/100.0),phi*phi*apowers
c     2               ,(2*ak*ak*phi*y(1)/grhom/omegam)**2*apowers
c                   enddo
c               endif

                       itau=itau+1
                     endif
                    endif
                    itf=itf+1
                    if (j.lt.nstep.and.itf.le.ntf.and.
     2                 atau0(j+1).gt.tautf(itf)) goto 101
                    endif
                   endif
                  end if
 10            end do
c     TESTING. Output averaged sources.
c               write(25,'(4E15.5)')ak0(ik),averaget,averagep,beta

               d(ik,nstep)=0.0d0
               dp(ik,nstep)=0.0d0
            end if
c
c     time loop for tensors if requested.
c
c     Tensors will only be calculated if k*tau< stsp.
c
            stpt=50.0d0
c
c
            if (itflag.ne.0) then
               call finitialt(yt,taustart0)
               tau=taustart0
               ind=1
c     Begin timestep loop.
c     dt contains the sources for the anisotropy and dtp
c     for the polarization. t means tensor.
               do 11 j=2,nstep
                  tauend=atau0(j)
                  if ((ak*tauend).gt.stpt) then
                     dt(ik,j)=0.0d0
                     dte(ik,j)=0.0d0
                     dtb(ik,j)=0.0d0
                  else
                     tol1=tol
                     call dverk(nvart,fderivst,tau,yt,tauend,tol1,
     &                    ind,c,nwt,wt)
                     call fderivst(nvart,tau,yt,ytprime)
                     call foutputt(nvart,yt,ytprime,j,tau0,tau,
     &                    dt(ik,j),dte(ik,j),dtb(ik,j))
                  end if
 11            end do

               dt(ik,nstep)=0.0d0
               dte(ik,nstep)=0.0d0
               dtb(ik,nstep)=0.0d0
            end if
c
c
c
 30      end do
c     TESTING: Close output of average
c         close(25)
      endif

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
#ifdef MAX

c     FOR PROJECT WITH MAX
!!      CHANGED BY MAX:
            nkt=int((log(akmaxt/(h0/100.))-log(1.d-5))*nlnkt)+1
            akdone=1.d-5
!! WANT akmaxt TO BE IN 1/Mpc, NOT h/Mpc.

#endif
         endif

         akgf=1.0d0
         if (akdone.gt.akgf) akgf=akdone*exp(1.5/nlnkt)

c     Loop over wavenumbers.
         do 31 ik=nk+1+myid,nkt,numprocs                !MPI modification
            ak=akdone*exp(dble(ik-nk)/nlnkt)
            if (lensflag.ne.0) itau=1

cxxxc     TIMING
cxxx      timeprev=actual
cxxx      actual=etime(tarray)
cxxx      write(50,*)actual-timeprev,'ik=',ik,ak
cxxx      write(*,*)actual-timeprev,'ik=',ik,ak

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
            taustart=0.001d0/ak
c     Make sure to start early in the radiation era.
            taustart=min(taustart,0.1d0)
c     Start when massive neutrinos are strongly relativistic.
            if (amnu.ne.0.0d0) then
               arel=1.d-3/amnu
               taustart=min(taustart,arel/adotrad)
            end if
            ak2=ak*ak
            call finitial(y,taustart)
            tau=taustart

#ifndef HISPDHIKMODES
            do 135 itf=1,ntf
              if ((tautf(itf)-tau).gt.1.0d-5) then
               call dverk(nvar,fderivs,tau,y,tautf(itf),tol,
     2              ind,c,nw,w)
              else
               tau=tautf(itf)
              endif
c
                     if(indtau(itf).ne.0) then
                      call outtransf(nvar,y,0.0d0,ik,indtau(itf))
                     else
                      if (lensflag.ne.0) then
                       call fderivs(nvar,tau,y,yprime)
                       call foutput(nvar,y,yprime,j,tau0,tau
     &                 ,temp,temp,temp,phi)
                       if (itau.eq.1) akl(ik)=ak
                       if (ik.eq.1) then 
                        atrg(itau,1)=tau
                        atrg(itau,2)=y(1)
                        atrg(itau,3)=tau
                       endif
                       phik(ik,itau)=phi
c                       write(17,*)y(1),tau,phi
c              if (itau.eq.ntau) then 
c                   do in=1,nn
c                     call powersflat(ak,in,apowers)
c                     write(21,'(5E14.5)')ak/(h0/100.0),phi*phi*apowers
c
c     2               ,(2*ak*ak*phi*y(1)/grhom/omegam)**2*apowers
c                   enddo
c              endif
                       itau=itau+1
                     endif
                    endif
 135        enddo
#else
            if(ak.le.akgf) then
               do 1352 itf=1,ntf
                if ((tautf(itf)-tau).gt.1.0d-5) then
                  call dverk(nvar,fderivs,tau,y,tautf(itf),tol,
     2              ind,c,nw,w)
                else
                 tau=tautf(itf)
                endif
                     if(indtau(itf).ne.0) then
                      call outtransf(nvar,y,0.0d0,ik,indtau(itf))
                     else
                      if (lensflag.ne.0) then
                       if (itf.eq.1) then 
                        call fderivs(nvar,tau,y,yprime)
                        call foutput(nvar,y,yprime,j,tau0,tau
     &                  ,temp,temp,temp,phi1)
                       endif
                       call fderivs(nvar,tau,y,yprime)
                       call foutput(nvar,y,yprime,j,tau0,tau
     &                 ,temp,temp,temp,phi)
                       if (itau.eq.1) akl(ik)=ak
                       if (ik.eq.1) then 
                        atrg(itau,1)=tau
                        atrg(itau,2)=y(1)
                        atrg(itau,3)=tau
                       endif
                       phik(ik,itau)=phi
c                       write(17,*)y(1),tau,phi
c              if (itau.eq.ntau) then 
c                   do in=1,nn
c                     call powersflat(ak,in,apowers)
c                     write(21,'(5E14.5)')ak/(h0/100.0),phi*phi*apowers
c     2               ,(2*ak*ak*phi*y(1)/grhom/omegam)**2*apowers
c                   enddo
c               endif
                       itau=itau+1
                     endif
                    endif
 1352          enddo

               closestknum=istore(1)-1
            else
               call dverk(nvar,fderivs,tau,y,tautf(1),tol,
     2              ind,c,nw,w)
                 if(indtau(itf).ne.0) then
                   call outtransf(nvar,y,0.0d0,ik,indtau(itf))
                 endif
               do 1353 itf=2,ntf
                 if(indtau(itf).ne.0) then
                  growthfactor=(omegac*power(6*closestknum+2,indtau(itf))
     $                 +omegab*power(6*closestknum+3,indtau(itf))
     $                 +omegan*power(6*closestknum+6,indtau(itf))
     $                 /(omegac*power(6*closestknum+2,1)
     $                 +omegab*power(6*closestknum+3,1)
     $                 +omegan*power(6*closestknum+6,1))
c     k
                  power(6*istore(indtau(itf))+1,itf)=
     1            power(6*istore(indtau(itf))+1,1)
c     all others are transfer functions
                  power(6*istore(itf)+2,indtau(itf))=
     $                 power(6*istore(indtau(itf))+2,1)*growthfactor
                  power(6*istore(indtau(itf))+3,indtau(itf))=
     $                 power(6*istore(indtau(itf))+3,1)*growthfactor
                  power(6*istore(indtau(itf))+4,indtau(itf))=
     $                 power(6*istore(indtau(itf))+4,1)*growthfactor
                  power(6*istore(indtau(itf))+5,indtau(itf))=
     $                 power(6*istore(indtau(itf))+5,1)*growthfactor
                  power(6*istore(indtau(itf))+6,indtau(itf))=
     $                 power(6*istore(indtau(itf))+6,1)*growthfactor
                  istore(indtau(itf)) = istore(indtau(itf)) +1
                else
                      if (lensflag.ne.0) then
                       call fderivs(nvar,tau,y,yprime)
                       call foutput(nvar,y,yprime,j,tau0,tau
     &                 ,temp,temp,temp,phihk)
                       if (itau.eq.1) akl(ik)=ak
                       if (ik.eq.1) then 
                        atrg(itau,1)=tau
                        atrg(itau,2)=y(1)
                        atrg(itau,3)=tau
                       endif
                       phik(ik,itau)=phihk*phi/phi1
c                       write(17,*)y(1),tau,phi
c              if (itau.eq.ntau) then 
c                   do in=1,nn
c                     call powersflat(ak,in,apowers)
c
c                     write(21,'(5E14.5)')ak/(h0/100.0),phi*phi*apowers
c     2               ,(2*ak*ak*phi*y(1)/grhom/omegam)**2*apowers
c
c                   enddo
c               endif
                       itau=itau+1
                     endif
                    endif

 1353          enddo
            endif
#endif
 31      end do
      endif

       if (lensflag.ne.0) ntf=ntf-ntau
#ifdef MPIBUILD
c
c       communicate values of sources to other processors
c
c
        if(myid.eq.0) write(*,*) 'starting communication'
c       write(*,*)myid,'has reached barrier'
      CALL MPI_BARRIER_(MPI_COMM_WORLD,IERR)
      CALL COMMUNICATE(d,dp,dkk,dt,dte,dtb,nstep,nk,itflag)
      CALL MPI_BARRIER_(MPI_COMM_WORLD,IERR)
#endif

c
c
      if(ict.ne.0) then
         call output_power(ntf,amnu)
#ifdef LYAOUTPUTS

#ifndef HISPDHIKMODES
         closestk=1.0e10
         closestknum=1

c     here we find k closest to k=akgf

         do itf = 0,istore(ntf)
            temp=abs(log(power(6*itf+1,ntf)/akgf*100.0/h0))
            if (temp.le.closestk) then
               closestknum=itf
               closestk=temp
            endif
         end do
#endif

         do itf=ntf,1,-1
            atf=1.0d0/(ztf(itf)+1.0d0)
            dtaudatf=dtauda(atf)
            htf=1.0d0/dtaudatf/atf**2*2.998d5
            if (ndyn.eq.0) then
               omegaq = omegav*atf**4*grhom*dtaudatf**2/3
            else
               omegaq=omegav*dynrho(atf)*atf**4*grhom*dtaudatf**2/3
            endif

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


c
cxxxc     TIMING
cxxx      timeprev=actual
cxxx      actual=etime(tarray)
cxxx      write(50,*)actual-timeprev,' END SOURCES, NOW Cls'
cxxx      write(*,*)actual-timeprev,' END SOURCES, NOW Cls'

c
c
c
c     if CMB calculations are requested, calculate the Cl by
c     integrating the sources over time and over k.
c
c
      if (ict.ne.1) then
c     get the interpolation matrix for the sources to interpolate them
c     for other k-values, scalar case.
         if (itflag.ne.2) then
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

#ifdef HP
         dk=1.1d0/tau0
         dk0=1.1d0/tau0
         no=700
         dlnk1=0.07d0
#else
         dk=2.5d0/tau0
         dk0=1.5d0/tau0
         no=700
         dlnk1=0.07d0
#endif
cxxx         write(*,*)'Enter accfact'
cxxx         read(*,*)accfact
cxxx         dk=dk/accfact
cxxx         dk0=dk0/accfact
cxxx         no=no*accfact
cxxx         dlnk1=dlnk1/accfact

         no1=int(log(10.0d0*dk0/akmin)/dlnk1)+1
         if (akmax.gt.(no*dk0)) then
            nko=int((akmax-no*dk0)/dk)+no
         else
            no=int((akmax-10.0d0*dk0)/dk0)+no1
            nko=no
         end if

         if (nko.gt.nkmax) then
            write(*,*)
     2         'Sorry, the arrays were dimensioned for a maximum of'
            write(*,*) nkmax, 'k modes.'
            write(*,*)'The model you requested needs',nko
            write(*,*)'Please make the arrays bigger by making '
            write(*,*)'nkmax bigger where it appears'
            stop
         end if

cxxxc     TIMING
cxxx        write(*,*)'ak integral',no1,no,nko
cxxx        write(50,*)'ak integral',no1,no,nko


         do 198 k=1,nko
            if (k.le.no) then
               if (k.le.no1) then
                  ak1(k)=10.0d0*dk0*exp(-(no1-k)*dlnk1)
                  dak1(k)=ak1(k)*dlnk1
               else
                  ak1(k)=ak1(no1)+(k-no1)*dk0
                  dak1(k)=dk0
               end if
            else
               ak1(k)=ak1(no)+(k-no)*dk
               dak1(k)=dk
            end if
 198     continue

c         do 199 k=2,(nko-1)
c 199     continue
         dak1(1)=0.5d0*dak1(1)
         dak1(no1)=0.5d0*(dak1(no1)+dk0)
         dak1(nko)=0.5d0*dak1(nko)
c
         klo=1
         khi=2

c     TESTING
c     To output average of the sources
c         open(unit=24,file='average_2.dat',status='unknown'
c     $        ,form='formatted')
c         rewind 24


         kkmin=1
         kkmax=nko

c     K split
         if ((kcutflag).eq.1) then
            do k=1,nko
               if (ak1(k).lt.(1.3d0*aksplit))kkmax=k+1
            end do
         end if
         if ((kcutflag).eq.-1) then
            do k=1,nko
               if (ak1(k).le.(0.3d0*aksplit))kkmin=k
            end do
         end if

c     TESTING
c         write(*,*)'kkmin,kkmax=',kkmin,kkmax

c     Begin k-loop
         do 200 k=kkmin+myid,kkmax,numprocs
            akt=ak1(k)

c     TESTING Output ks where the source is calculated.
c     write(445,*)k,akt

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

            nstps=0
            nstpt=0

            if (itflag.ne.2) then
               do 250 j=1,l0
                  dl(j)=0.0d0
                  dl2(j)=0.0d0
                  dl3(j)=0.0d0
                  dpl(j)=0.0d0
                  dpl2(j)=0.0d0
                  dpl3(j)=0.0d0
                  dkl(j)=0.0d0
                  dkl2(j)=0.0d0
                  dkl3(j)=0.0d0
 250           continue
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

               s2(1)=0.0d0
               sp2(1)=0.0d0
               sk2(1)=0.0d0
               do 304 i=2,nstps
                  s2(i)=a0*d(klo,i)+b0*d(khi,i)+((a0**3-a0)
     $                 *dpr(klo,i)
     &                 +(b0**3-b0)*dpr(khi,i))*ho*ho/6.d0
                  sp2(i)=a0*dp(klo,i)+b0*dp(khi,i)+((a0**3-a0)
     &                 *dppr(klo,i)+(b0**3-b0)*dppr(khi,i))
     $                 *ho*ho/6.0d0
                  sk2(i)=a0*dkk(klo,i)+b0*dkk(khi,i)+((a0**3-a0)
     &                 *dkpr(klo,i)+(b0**3-b0)*dkpr(khi,i))
     $                 *ho*ho/6.0d0

c     TESTING. Compute averaged sources for each k in the time
c     integral to compare with those calculated earlier.
c                  averaget=averaget+s2(i)*dtau1(i)
c                  averagep=averagep+sp2(i)*dtau1(i)

 304           continue

c     TESTING. Output averaged sources.
c               write(24,'(4E15.5)')akt,averaget,averagep,beta

               s2(nstps+1)=0.0d0
               sp2(nstps+1)=0.0d0
               sk2(nstps+1)=0.0d0
            end if

c     If tensors wanted
            if (itflag.ne.0) then
               do 260 j=1,l0
                  dtl(j)=0.0d0
                  dtl2(j)=0.0d0
                  dtl3(j)=0.0d0
                  dtel2(j)=0.0d0
                  dtbl2(j)=0.0d0
 260           continue
c
c     Interpolating the tensor source as a function of time for the present
c     wavelength.
               st2(1)=0.0d0
               ste2(1)=0.0d0
               stb2(1)=0.0d0
               nstpt=2
               do 306 i=2,nstep
                  xf=akt*(tau0-atau0(i))
                  if (((akt*atau0(i)).lt.stpt).
     $                 and.(xf.gt.1.0d-8)) then
                     nstpt=i
                     st2(i)=a0*dt(klo,i)+b0*dt(khi,i)+((a0**3-a0)
     &                    *dtpr(klo,i)+(b0**3-b0)*dtpr(khi,i))
     $                    *ho*ho/6.0d0
                     ste2(i)=a0*dte(klo,i)+b0*dte(khi,i)+((a0**3-a0)
     &                    *dtepr(klo,i)+(b0**3-b0)*dtepr(khi,i))
     $                    *ho*ho/6.0d0
                     stb2(i)=a0*dtb(klo,i)+b0*dtb(khi,i)+((a0**3-a0)
     &                    *dtbpr(klo,i)+(b0**3-b0)*dtbpr(khi,i))
     $                    *ho*ho/6.0d0
                  else
                     st2(i)=0.0d0
                     ste2(i)=0.0d0
                     stb2(i)=0.0d0
                  end if
 306           end do
               nstpt=max(nstpt,n1)
               st2(nstpt+1)=0.0d0
               ste2(nstpt+1)=0.0d0
               stb2(nstpt+1)=0.0d0
            end if

c     Findind the position in the xx table for the x correponding to each
c     timestep
            do 405 i=1,max(nstps,nstpt)+2
               xf=abs(akt*(tau0-atau0(i)))
               if (xf.le.25.0) then
                  if (xf.le.5.0) then
                     de2=10.0d0*xf+1.0d0
                  else
                     de2=(xf-5.0d0)*5.0d0+51.0d0
                  end if
               else
                  de2=(xf-25.0d0)+151.0d0
               end if
               mxx(i)=int(de2)
               mxx(i)=max(mxx(i),1)
 405        end do

c
c     Begin l and  time-loop to integrate scalar perturbations.
c
c     Determining ranges of integration

            if (itflag.ne.2) then
               do 310 j=1,l0
                  xlim=0.05d0*l(j)
                  xlim=max(xlim,xlimmin)
                  xlim=l(j)-xlim
                  xlmax1=80.0d0*l(j)
                  xlmax2=dble(min(2*l(j),kmax))
c     TESTING
c                  xlmax2=2.0d0*xlmax2

                  tmin=tau0-xlmax1/akt
                  tmax=tau0-xlim/akt
                  tmax=min(tau0,tmax)
                  tmin=max(atau0(2),tmin)
c     TESTING
c                  tmin=atau0(2)
c
                  if (tmax.lt.atau0(2)) goto 400

                  if (zri.eq.0.0) then

                     if (tmin.lt.taurend) then
                        nstart1=2
                     else
                        nstart1=n1+int(log(tmin/taurend)/dlntau0)
                     end if
                     if (tmax.lt.taurend) then
                        nstop1=n1
                        nstop1a=n1
                     else
                        nstop1=n1+int(log(tmax/taurend)/dlntau0)
                        nstop1=min(nstop1,nstps)
                        if ((akt*dtau1(nstop1)).gt.1.0) then
                           nstop1a=max(n1,n1-int(log(akt*dlntau0
     &                          *taurend)/dlntau0))
                        else
                           nstop1a=nstop1
                        end if
                     end if

                  else

                     if (tmin.lt.taurend) then
                        nstart1=2
                        nstart2=j2ri1
                     else
                        if (tmin.lt.atau0(nriend)) then
                           nstart1=min(n1+int(log(tmin/taurend)
     &                          /dlntau0),j2ri1)
                           nstart2=j2ri1
                        else
                           nstart1=nstep+1
                           nstart2=nriend+int(log(tmin/atau0(nriend)
     &                          /dlntau0))
                        end if
                     end if

                     if (tmax.lt.taurend) then
                        nstop1=n1
                        nstop2=0
                     else
                        if (tmax.lt.atau0(j2ri1)) then
                           nstop1=n1+int(log(tmax/taurend)/dlntau0)
                           nstop2=0
                        else
                           nstop1=j2ri1-1
                           nstop2=max(nriend,nriend+int(log(tmax
     &                          /atau0(nriend))/dlntau0))

                           nstop2=min(nstop2,nstps)
                           if ((akt*dtau1(nstop2)).gt.1) then
                              nstop2a=max(nriend,nriend-int(log(akt
     &                             *dlntau0*atau0(nriend))/dlntau0))
                           else
                              nstop2a=nstop2
                           end if

                        end if

                        nstop1=min(nstop1,j2ri1-1)
                        if ((akt*dtau1(nstop1)).gt.1) then
                           nstop1a=max(n1,n1-int(log(akt
     &                          *dlntau0*taurend)/dlntau0))
                        else
                           nstop1a=nstop1
                        end if

                     end if
                  end if


c     Integration before reionization.
c
c     Interpolating jls at points where the
c     sources are recorded.
                  do i=nstart1,nstop1a
                     xf=akt*(tau0-atau0(i))
                     m2=mxx(i)
                     h2=xx(m2+1)-xx(m2)
                     a2=(xx(m2+1)-xf)/h2
                     b2=(xf-xx(m2))/h2
                     ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)+((a2**3-a2)
     &                    *ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))
     $                    *(h2*h2)/6.d0

                     dl2(j)=dl2(j)+s2(i)*ajl0*dtau2(i)
                     dpl2(j)=dpl2(j)+sp2(i)*ajl0*dtau2(i)
                     dkl2(j)=dkl2(j)+sk2(i)*ajl0*dtau2(i)
                  end do

c     Interpolating sources at points where the
c     jls are recorded.

                  do i=nstop1a+1,nstop1
                     xf=akt*(tau0-atau0(i))
                     m2=mxx(i)
                     dtau3=dtau1(i)*abs(s2(i)
     $                    /(s2(i+1)-s2(i)+1.0d-10))
c     TESTING
c                     dtau3=dtau3/20.0d0

                     if ((xf.lt.xlmax2).or.
     $                    ((akt*dtau3).lt.1.0)) then
                        xi=xf-akt*dtau1(i)
                        m1=mxx(i+1)
                        ddt1=s2(i)
                        ddt2=s2(i+1)
                        ddp1=sp2(i)
                        ddp2=sp2(i+1)
                        ddk1=sk2(i)
                        ddk2=sk2(i+1)
                        do lx=m1+1,m2
                           x=xx(lx)
                           ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2
                           ddp=(ddp1-ddp2)*(x-xi)/(xf-xi)+ddp2
                           ddk=(ddk1-ddk2)*(x-xi)/(xf-xi)+ddk2
                           dl3(j)=dl3(j)+ajl(lx,j)*ddt*dxx(lx)
                           dpl3(j)=dpl3(j)+ajl(lx,j)*ddp*dxx(lx)
                           dkl3(j)=dkl3(j)+ajl(lx,j)*ddk*dxx(lx)
                        end do
                     end if
                  end do

c     Integration after reionization
c     Interpolating jls at points where the
c     sources are recorded.

                  if (zri.ne.0.0) then
                     do i=nstart2,nstop2a
                        xf=akt*(tau0-atau0(i))
                        m2=mxx(i)
                        h2=xx(m2+1)-xx(m2)
                        a2=(xx(m2+1)-xf)/h2
                        b2=(xf-xx(m2))/h2
                        ajl0=a2*ajl(m2,j)
     $                       +b2*ajl(m2+1,j)+((a2**3-a2)*
     &                       ajlpr(m2,j)+(b2**3-b2)*ajlpr(m2+1,j))
     &                       *(h2*h2)/6.d0
                        dl2(j)=dl2(j)+s2(i)*ajl0*dtau2(i)
                        dpl2(j)=dpl2(j)+sp2(i)*ajl0*dtau2(i)
                        dkl2(j)=dkl2(j)+sk2(i)*ajl0*dtau2(i)
                     end do

c     Interpolating sources at points where the
c     jls are recorded.

                     do i=nstop2a+1,nstop2
                        xi=xf-akt*dtau1(i)
                        xf=akt*(tau0-atau0(i))
                        m2=mxx(i)
                        dtau3=dtau1(i)*abs(s2(i)
     $                       /(s2(i+1)-s2(i)+1.0d-10))
c     TESTING
c                     dtau3=dtau3/2.0d0

                        if ((xf.lt.xlmax2).or
     $                       .((akt*dtau3).lt.1)) then
                           m1=mxx(i+1)
                           ddt1=s2(i)
                           ddt2=s2(i+1)
                           ddp1=sp2(i)
                           ddp2=sp2(i+1)
                           ddk1=sk2(i)
                           ddk2=sk2(i+1)
                           do lx=m1+1,m2
                              x=xx(lx)
                              ddt=(ddt1-ddt2)*(x-xi)/(xf-xi)+ddt2
                              ddp=(ddp1-ddp2)*(x-xi)/(xf-xi)+ddp2
                              ddk=(ddk1-ddk2)*(x-xi)/(xf-xi)+ddk2
                              dl3(j)=dl3(j)+ajl(lx,j)*ddt*dxx(lx)
                              dpl3(j)=dpl3(j)+ajl(lx,j)*ddp*dxx(lx)
                              dkl3(j)=dkl3(j)+ajl(lx,j)*ddk*dxx(lx)
                           end do
                        end if
                     end do
                  end if

                  dl(j)=dl2(j)+dl3(j)/akt
                  dpl(j)=dpl2(j)+dpl3(j)/akt
                  dkl(j)=dkl2(j)+dkl3(j)/akt
 310           continue

 400           continue
c     TESTING. Output the result of line of sight integral
c     for several ls as a function of k.
c     write(350,'(5E15.5)')akt,dl(1),dl(2),dl(3),dl(4)
c     write(351,'(5E15.5)')akt,dl(5),dl(6),dl(7),dl(8)
c     write(352,'(5E15.5)')akt,dl(9),dl(10),dl(11),dl(12)
c     write(353,'(5E15.5)')akt,dl(13),dl(14),dl(15),dl(16)
c     write(354,'(5E15.5)')akt,dl(17),dl(18),dl(19),dl(20)
c     write(355,'(5E15.5)')akt,dpl(1),dpl(2),dpl(3),dpl(4)
c     write(356,'(5E15.5)')akt,dpl(5),dpl(6),dpl(7),dpl(8)
c     write(357,'(5E15.5)')akt,dpl(9),dpl(10),dpl(11),dpl(12)
c     write(358,'(5E15.5)')akt,dpl(13),dpl(14),dpl(15),dpl(16)
c     write(359,'(5E15.5)')akt,dpl(17),dpl(18),dpl(19),dpl(20)

c     Adding to calculate the integral over k.
c     Scalar case
c
               do in=1,nn
                  do j=1,l0
                     call powersflat(akt,in,apowers)
                     apowers=apowers/akt
                     ckj=apowers*dl(j)*dl(j)*dak1(k)
                     cl(j,in)=cl(j,in)+ckj
                     cpkj=apowers*dpl(j)*dpl(j)*dak1(k)
                     cpl(j,in)=cpl(j,in)+cpkj
                     cckj=apowers*dl(j)*dpl(j)*dak1(k)
                     ccl(j,in)=ccl(j,in)+cckj
                     ckkkj=apowers*dkl(j)*dkl(j)*dak1(k)
                     ckkl(j,in)=ckkl(j,in)+ckkkj
                     ctkkj=apowers*dl(j)*dkl(j)*dak1(k)
                     ctkl(j,in)=ctkl(j,in)+ctkkj

c     TESTING. Keep integrand of cls for each k.
c     ctt(j)=ckj
c     cpp(j)=cpkj

                  end do
               end do

c     TESTING. Output Integrand of cl as a function of k for
c     several ls.
c     write(450,'(5E15.5)')akt,ctt(1),ctt(2),ctt(3),ctt(4)
c     write(451,'(5E15.5)')akt,ctt(5),ctt(6),ctt(7),ctt(8)
c     write(452,'(5E15.5)')akt,ctt(9),ctt(10),ctt(11),ctt(12)
c     write(453,'(5E15.5)')akt,ctt(13),ctt(14)
c     $            ,ctt(15),ctt(16)
c     write(454,'(5E15.5)')akt,ctt(17),ctt(18)
c     $            ,ctt(19),ctt(20)
c     write(455,'(5E15.5)')akt,cpp(1),cpp(2),cpp(3),cpp(4)
c     write(456,'(5E15.5)')akt,cpp(5),cpp(6),cpp(7),cpp(8)
c     write(457,'(5E15.5)')akt,cpp(9),cpp(10),cpp(11),cpp(12)
c     write(458,'(5E15.5)')akt,cpp(13),cpp(14)
c     $            ,cpp(15),cpp(16)
c     write(459,'(5E15.5)')akt,cpp(17),cpp(18)
c     $            ,cpp(19),cpp(20)
c

            end if

c
c     Begin l and  time-loop to integrate tensor perturbations.
c
c     Finding the ranges of integration.

            if (itflag.ne.0) then
               do 530 j=1,l0
                  xlim=0.05d0*l(j)
                  xlim=max(xlim,xlimmin)
                  xlim=l(j)-xlim
                  xlmax1=80.0d0*l(j)
c     xlmax2=dble(min(2*l(j),kmax))
                  tmin=tau0-xlmax1/akt
                  tmax=tau0-xlim/akt
                  tmax=min(tau0,tmax)
                  tmin=max(atau0(2),tmin)

c
                  if (tmax.lt.atau0(2)) goto 590

                  if (zri.eq.0.0) then

                     if (tmin.lt.taurend) then
                        nstart1=2
                     else
                        nstart1=n1+int(log(tmin/taurend)/dlntau0)
                     end if
                     if (tmax.lt.taurend) then
                        nstop1=n1
                        nstop1a=n1
                     else
                        nstop1=n1+int(log(tmax/taurend)/dlntau0)
                        nstop1=min(nstop1,nstpt)
                        if ((akt*dtau1(nstop1)).gt.1) then
                           nstop1a=max(n1,n1-int(log(akt
     &                          *dlntau0*taurend)/dlntau0))
                        else
                           nstop1a=nstop1
                        end if
                     end if

                  else

                     if (tmin.lt.taurend) then
                        nstart1=2
                        nstart2=j2ri1
                     else
                        if (tmin.lt.atau0(nriend)) then
                           nstart1=min(n1+int(log(tmin/taurend)
     &                          /dlntau0),j2ri1)
                           nstart2=j2ri1
                        else
                           nstart1=nstep+1
                           nstart2=nriend+int(log(tmin/atau0(nriend)
     &                          /dlntau0))
                        end if
                     end if

                     if (tmax.lt.taurend) then
                        nstop1=n1
                        nstop2=0
                     else
                        if (tmax.lt.atau0(j2ri1)) then
                           nstop1=n1+int(log(tmax/taurend)/dlntau0)
                           nstop2=0
                        else
                           nstop1=j2ri1-1
                           nstop2=max(nriend,nriend+int(log(tmax
     &                          /atau0(nriend))/dlntau0))

                           nstop2=min(nstop2,nstpt)
                           if ((akt*dtau1(nstop2)).gt.1) then
                              nstop2a=max(nriend,nriend-int(log(akt
     &                             *dlntau0*atau0(nriend))/dlntau0))
                           else
                              nstop2a=nstop2
                           end if

                        end if

                        nstop1=min(nstop1,j2ri1-1)
                        if ((akt*dtau1(nstop1)).gt.1) then
                           nstop1a=max(n1,n1-int(log(akt
     &                          *dlntau0*taurend)/dlntau0))
                        else
                           nstop1a=nstop1
                        end if

                     end if
                  end if

c     Integration before reionization.
c
c     Interpolating jls at points where the
c     sources are recorded.

                  do i=nstart1,nstop1a
                     xf=akt*(tau0-atau0(i))
                     if ((j.eq.1).and.(xf.lt.1.0)
     $                    .and.(xf.gt.1.0d-7)) then
                        ajl0=(3.d0*(sin(xf)/xf-cos(xf))/xf
     $                       - sin(xf))/xf
                     else
                        m2=mxx(i)
                        h2=xx(m2+1)-xx(m2)
                        a2=(xx(m2+1)-xf)/h2
                        b2=(xf-xx(m2))/h2
                        ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)
     $                       +((a2**3-a2)
     &                       *ajlpr(m2,j)+(b2**3-b2)
     $                       *ajlpr(m2+1,j))
     &                       *(h2*h2)/6.d0
                    endif
                    dtl2(j)=dtl2(j)+st2(i)*ajl0*dtau2(i)
                    dtel2(j)=dtel2(j)+ste2(i)*ajl0*dtau2(i)
                    dtbl2(j)=dtbl2(j)+stb2(i)*ajl0*dtau2(i)
                 end do

c     Interpolating sources at points where the
c     jls are recorded.

                  do i=nstop1a+1,nstop1
                     xf=akt*(tau0-atau0(i))
                     m2=mxx(i)
                     xi=xf-akt*dtau1(i)
                     m1=mxx(i+1)
                     dtdt1=st2(i)
                     dtdt2=st2(i+1)
                     dtde1=ste2(i)
                     dtde2=ste2(i+1)
                     dtdb1=stb2(i)
                     dtdb2=stb2(i+1)
                     do lx=m1+2,m2+1
                        x=xx(lx)
                        dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
                        dtde=(dtde1-dtde2)*(x-xi)/(xf-xi)+dtde2
                        dtdb=(dtdb1-dtdb2)*(x-xi)/(xf-xi)+dtdb2
                        dtl3(j)=dtl3(j)+ajl(lx,j)*dtdt*dxx(lx)
                        dtel3(j)=dtel3(j)+ajl(lx,j)*dtde*dxx(lx)
                        dtbl3(j)=dtbl3(j)+ajl(lx,j)*dtdb*dxx(lx)
                     end do
                  end do

c     Integration after reionization
                  if (zri.ne.0.0) then
c     Interpolating jls at points where the
c     sources are recorded.
c
                     do i=nstart2,nstop2a
                        xf=akt*(tau0-atau0(i))

                        if ((j.eq.1).and.(xf.lt.1.0)
     $                       .and.(xf.gt.1.0d-7)) then
                           ajl0=(3.d0*(sin(xf)/xf-cos(xf))/xf
     $                          - sin(xf))/xf
                        else
                           m2=mxx(i)
                           h2=xx(m2+1)-xx(m2)
                           a2=(xx(m2+1)-xf)/h2
                           b2=(xf-xx(m2))/h2
                           ajl0=a2*ajl(m2,j)+b2*ajl(m2+1,j)
     $                          +((a2**3-a2)
     &                          *ajlpr(m2,j)+(b2**3-b2)
     $                          *ajlpr(m2+1,j))
     &                          *(h2*h2)/6.d0
                        endif

                        dtl2(j)=dtl2(j)+st2(i)*ajl0*dtau2(i)
                        dtel2(j)=dtel2(j)+ste2(i)*ajl0*dtau2(i)
                        dtbl2(j)=dtbl2(j)+stb2(i)*ajl0*dtau2(i)
                     end do

c     Interpolating sources at points where the
c     jls are recorded.

                     do i=nstop2a+1,nstop2
                        xi=xf-akt*dtau1(i)
                        xf=akt*(tau0-atau0(i))
                        m2=mxx(i)
                        m1=mxx(i+1)
                        dtdt1=st2(i)
                        dtdt2=st2(i+1)
                        dtde1=ste2(i)
                        dtde2=ste2(i+1)
                        dtdb1=stb2(i)
                        dtdb2=stb2(i+1)
                     do lx=m1+2,m2+1
                        x=xx(lx)
                        dtdt=(dtdt1-dtdt2)*(x-xi)/(xf-xi)+dtdt2
                        dtde=(dtde1-dtde2)*(x-xi)/(xf-xi)+dtde2
                        dtdb=(dtdb1-dtdb2)*(x-xi)/(xf-xi)+dtdb2
                        dtl3(j)=dtl3(j)+ajl(lx,j)*dtdt*dxx(lx)
                        dtel3(j)=dtel3(j)+ajl(lx,j)*dtde*dxx(lx)
                        dtbl3(j)=dtbl3(j)+ajl(lx,j)*dtdb*dxx(lx)
                     end do
                   end do

                  end if

                  dtl(j)=dtl2(j)+dtl3(j)/akt
                  dtel(j)=dtel2(j)+dtel3(j)/akt
                  dtbl(j)=dtbl2(j)+dtbl3(j)/akt
 530           continue

 590           continue


c     Tensor case.
c
c     Adding to calculate the integral over k.
c     Tensor case

               do in=1,nn
                  do j=1,l0
                     call powertflat(akt,in,apowert)
                     apowert=apowert/akt


                     ctkj=apowert*dtl(j)*dtl(j)*dak1(k)
                     ctl(j,in)=ctl(j,in)+ctkj

                     ctekj=apowert*dtel(j)*dtel(j)*dak1(k)
                     ctel(j,in)=ctel(j,in)+ctekj

                     ctbkj=apowert*dtbl(j)*dtbl(j)*dak1(k)
                     ctbl(j,in)=ctbl(j,in)+ctbkj

                     ctckj=apowert*dtl(j)*dtel(j)*dak1(k)
                     ctcl(j,in)=ctcl(j,in)+ctckj

                  end do
               end do
c
            end if

 200     continue

#ifdef MPIBUILD
      CALL TRANSFER_CL(cl,cpl,ccl,ckkl,ctkl,ctl,ctel,ctbl,ctcl,
     >                  itflag,nn,l0)
      IF(MYID.ne.0) return
#endif

c     TESTING: Close file for output average
c         close(24)

cxxxc     TIMING
cxxx      timeprev=actual
cxxx      actual=etime(tarray)
cxxx      write(50,*)actual-timeprev,' END Cls, NOW FINAL'
cxxx      write(*,*)actual-timeprev,' END Cls, NOW FINAL'

c
c
c
c     Final calculations for CMB output.
c
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
                  ccl(j,in)=2.0d0*sqrt(ctnorm)
     &                 *ccl(j,in)*dble(l(j)*(l(j)+1))
                  ckkl(j,in)=2.0d0
     &                 *ckkl(j,in)*dble(l(j)*(l(j)+1))
                  ctkl(j,in)=2.0d0
     &                 *ctkl(j,in)*dble(l(j)*(l(j)+1))

c     TESTING. Output cls that are really calculated, from which
c     interpolation table is constructed.
c     write(11,*)l(j),cl(j,in),cpl(j,in)

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
c     Tensor Case
         if (itflag.ne.0) then
c     Normalization
            do 610 j=1,l0
               ctnorm=dble(l(j)*(l(j)-1.0d0)*(l(j)+1)*(l(j)+2))
               do 640 in=1,nn
                  ctl(j,in)=ctnorm*ctl(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctel(j,in)=ctel(j,in)
     $                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctbl(j,in)=ctbl(j,in)
     $                 *dble(l(j)*(l(j)+1))/8.0d0
                  ctcl(j,in)=sqrt(ctnorm)*ctcl(j,in)
     &                 *dble(l(j)*(l(j)+1))/8.0d0

c     TESTING. Output cls that are really calculated, from which
c     interpolation table is constructed.
c                  write(11,*)l(j),ctl(j,in),ctel(j,in)

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
                  clint=a0*cl(llo,in)+b0*cl(lhi,in)+((a0**3-a0)*
     &                 clpr(llo,in)+(b0**3-b0)
     $                 *clpr(lhi,in))*ho*ho /6.d0
                  cplint=a0*cpl(llo,in)+b0*cpl(lhi,in)+((a0**3-a0)
     &                 *cplpr(llo,in)+(b0**3-b0)
     $                 *cplpr(lhi,in))*ho*ho /6.d0
                  cclint=a0*ccl(llo,in)+b0*ccl(lhi,in)+((a0**3-a0)
     &                 *cclpr(llo,in)+(b0**3-b0)
     $                 *cclpr(lhi,in))*ho*ho /6.d0
                  ckklint=a0*ckkl(llo,in)+b0*ckkl(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ckklpr(llo,in)+(b0**3-b0)
     $                 *ckklpr(lhi,in))*ho*ho /6.d0
                  ctklint=a0*ctkl(llo,in)+b0*ctkl(lhi,in)+
     $                 ((a0**3-a0)
     &                 *ctklpr(llo,in)+(b0**3-b0)
     $                 *ctklpr(lhi,in))*ho*ho /6.d0
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

cxxxc     TIMING
cxxx      timeprev=actual
cxxx      actual=etime(tarray)
cxxx      write(50,*)actual-timeprev,' END FINAL, NOW BACK'
cxxx      write(*,*)actual-timeprev,' END FINAL, NOW BACK'

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fderivs(n,x,y,yprime)

c  Evaluate the time derivatives of the perturbations.
c
        implicit double precision (a-h,o-z)
        integer initfl

        integer ndyn
        double precision wdyn

        dimension y(n),yprime(n)
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
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
        common /store/ denl,dlfdlq
c
        common /epsilon/ epsw
c
        common /qparm/ wdyn,ndyn

c ep is used to stop the tight coupling approximation.
        if (ak.gt.epsw) then
           ep=ep0
        else
           ep=0.5d0*ep0
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
     $       *annur+annunr*rhonu*grhonr)/a2
     2      +grhom*omegavdyn*a2

#ifdef DIM

c energy density for  5 dimensions
        grho=(sqrt(grho)+sqrt(omegav*grhom*a2))**2

#endif

        adotoa=sqrt(grho/3.0d0)
        yprime(1)=adotoa*a
        gpres=((grhog+grhor*annur)/3.0d0+grhonr*annunr*pnu)/a2
     2      +weos*grhom*omegavdyn*a2

        if (ndyn.eq.1.or.ndyn.eq.2) then

c dark energy perturbations
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

        etadot=0.5d0*dgtheta/ak2
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

        rgrho=sqrt(grhom*omegav*grhom*omegav
        1        +grhom*(omegab+omegac)/(a2*a)+
     2  (grhog+grhor*(annur+annunr*rhonu))/(a2*a2))



          adotdota1=sgrhooa2/rgrho*(grhom*(omegab+omegac)*(-3/(2*a))
     2          +(grhog+grhor*(annur+annunr*rhonu))*(-2/a2)
     3          +grhor*annunr*rhonudot/(2*a2*adotoa))
     4          +2*grho

        adotdota=adotdota1/3.0d0
#else

        adotdota=0.5d0*(adotoa*adotoa-gpres)

#endif

c  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
          slip=(2.0d0*pb43/(1.0d0+pb43)+dlxedla)
     $       *adotoa*(thetab-thetag)
     2     +1.0d0/opac*(-adotdota*thetab-adotoa*ak2*0.5d0*deltag
     3     +ak2*(cs2*deltabdot-0.25d0*deltagdot))/(1.0d0+pb43)
c  First-order approximation to baryon velocity.
          thetabdot=thetabdot+pb43/(1.0d0+pb43)*slip
        end if
        yprime(7)=thetabdot
c  Photon total intensity and polarization equations of motion.
        yprime(8)=deltagdot
        thetagdot=(-thetabdot-adotoa*thetab+ak2*cs2*deltab)/pb43
     2   +ak2*(0.25d0*deltag-shearg)
        yprime(9)=thetagdot
        if (tcpa.eq.1) then
c  Treat baryons and photons as uncoupled.
          yprime(10)=8.0d0/15.0d0*thetag-0.6d0*ak*y(11)-opac*y(10)
     2        +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot
     $          +0.1d0*opac*polter
c  Polarization equations for l = 0, 1, 2.
          yprime(9+lmax)=-ak*y(10+lmax)-opac*y(9+lmax)
     $         +0.5d0*opac*polter
          yprime(10+lmax)=ak/3.0d0*(y(9+lmax)-2.0d0*y(11+lmax))
     2             -opac*y(10+lmax)
          yprime(11+lmax)=ak*(0.4d0*y(10+lmax)-0.6d0*y(12+lmax))
     2             -opac*y(11+lmax)+0.1d0*opac*polter
            do 10 l=3,lmax-1
            yprime(8+l)=ak*denl(l)*(l*y(7+l)-(l+1)*y(9+l))
     $              -opac*y(8+l)
            yprime(9+lmax+l)=ak*denl(l)*(l*y(8+lmax+l)-(l+1)*
     2                 y(10+lmax+l))-opac*y(9+lmax+l)
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
c       yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)
c     -opac*y(9+2*lmax)
        yprime(8+lmax)=ak*y(7+lmax)-(lmax+1)/tau*y(8+lmax)
     2                -opac*y(8+lmax)
        yprime(9+2*lmax)=ak*y(8+2*lmax)-(lmax+1)/tau*y(9+2*lmax)
     2                -opac*y(9+2*lmax)
c  Massless neutrino equations of motion.
        deltardot=4.0d0/3.0d0*(-thetar-0.5d0*hdot)
        yprime(10+2*lmax)=deltardot
        thetardot=ak2*(0.25d0*deltar-shearr)
        yprime(11+2*lmax)=thetardot
        yprime(12+2*lmax)=8.0d0/15.0d0*thetar-0.6d0*ak*y(13+2*lmax)
     2             +4.0d0/15.0d0*hdot+8.0d0/5.0d0*etadot
          do 20 l=3,lmaxnr-1
          yprime(10+2*lmax+l)=ak*denl(l)*(l*y(9+2*lmax+l)
     2                    -(l+1)*y(11+2*lmax+l))
20      continue
c  Truncate moment expansion.
c       yprime(10+2*lmax+lmaxnr)=ak*lmax*
c     y(9+2*lmax+lmaxnr)/(2*lmaxnr+1)
        yprime(10+2*lmax+lmaxnr)=ak*y(9+2*lmax+lmaxnr)
     $       -(lmaxnr+1)/tau
     7  *y(10+2*lmax+lmaxnr)
c

       if(ndyn.eq.1.or.ndyn.eq.2) then

        call dyn_phi(a,hdot,grho,gpres,rc_phi,rc_psi,
     2  rc_dphi,rc_dpsi)         
        yprime(n-1) = rc_dphi
        yprime(n)   = rc_dpsi

       endif

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
          yprime(ind)=akv(i)*(y(ind-nqmax)-2*y(ind+nqmax))/3
          ind=iq2+i-1
          yprime(ind)=akv(i)*(2*y(ind-nqmax)-3*y(ind+nqmax))/5
     2           -(hdot/15.0d0+2.0d0/5.0d0*etadot)*dlfdlq(i)
          ind=10+2*lmax+lmaxnr+i+lmaxnu*nqmax
c  Truncate moment expansion.
c         yprime(ind)=akv*lmaxnu*y(ind-nqmax)/(2*lmaxnu+1)
          yprime(ind)=akv(i)*y(ind-nqmax)-(lmaxnu+1)/tau*y(ind)
30      continue
          do 50 l=3,lmaxnu-1
            do 40 i=1,nqmax
            ind=10+2*lmax+lmaxnr+i+l*nqmax
          yprime(ind)=akv(i)*denl(l)*
     $           (l*y(ind-nqmax)-(l+1)*y(ind+nqmax))
40        continue
50      continue

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine foutput(n,y,yprime,j,tau0,tau,d,dp,dk,phi)
        implicit double precision (a-h,o-z)
        integer n,j
c
        include 'cmbfast.inc'

c       parameter (nstep0=2400)
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     $       grhonr
        common /out1/ adotoa,hdot,dgshear,rhonu,shearnu
c
        dimension y(n),yprime(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        save /visib/
        common /maxreion/ taurmax,armax,jrmax

c
        x=ak*(tau0-tau)
        a=y(1)
        a2=a*a
c       ahdot=y(2)
        eta=y(3)
        etadot=yprime(3)
        alpha=(hdot+6*etadot)/(2.0d0*ak2)
        alphadot=-3*dgshear/(2.0d0*ak2)+eta-2.0d0*adotoa*alpha
c  Baryons.
        deltab=y(6)
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
c        coupl=8.0d0*(thetag+ak2*alpha)/15.0d0-ak*0.6d0
c     2      *(y(11)+y(10+lmaxg)+y(12+lmaxg))
        coupldot=8.0d0*(thetagdot+ak2*alphadot)/15.0d0
     2      -ak*0.6d0*(yprime(11)+yprime(10+lmaxg)+yprime(12+lmaxg))
        polterdot=yprime(10)+yprime(9+lmaxg)+yprime(11+lmaxg)
        polterddot=coupldot-0.3d0*(dopac(j)*polter+opac(j)*
     $       polterdot)
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
        grhodot=grhodot-((1.0d0+3.0d0*weos)*grhom*
     $       omegavdyn*a2)*adotoa

#ifdef DIM

        grho=(sqrt(grhom*(omegac+omegab)/a+(grhog+grhor
     1          *(annur+annunr*rhonu))/a2
     2      +grhom*omegav*a2)
     3      +sqrt(omegav*grhom)*a)**2

                sgrhooa2=sqrt(grho/a2)

        rgrho=sqrt(grhom*omegav+grhom*(omegab+omegac)/(a2*a)+
     2  (grhog+grhor*(annur+annunr*rhonu))/(a2*a2))


        grhodot=sgrhooa2/rgrho*(grhom*(omegab+omegac)*(-3/a)+
     2          (grhog+grhor*(annur+annunr*rhonu))*(-4/a2)+
     3          grhor*annunr*rhonudot/(a2*adotoa))+
     4          2*grho

#endif

        adotoadot=grhodot/(6*adotoa)
c  Derivative of the shear
        dgsheardot=4.0d0/3.0d0*(grhog*sheargdot+annur*grhor
     2  *shearrdot)/a2-2.0d0*adotoa*dgshear+
     $       annunr*grhonr*shearnudot/a2
c  Calculation of the sources
        alphaddot=-3*dgsheardot/(2.0d0*ak2)+etadot
     2      -2.0d0*adotoadot*alpha-2.0d0*adotoa*alphadot
c
        s1=etadot+alphaddot
        s2=2*alphadot
c
c     TESTING
c        s1=0.0d0

        d=expmmu(j)*s1+vis(j)*(0.25d0*deltag
     2     +s2+polter/16+thetabdot/ak2+3.0d0/16.0d0
     3     *polterddot/ak2)+dvis(j)*(alpha+thetab/ak2+3.0d0/8.0d0
     4     *polterdot/ak2)+ddvis(j)*3.0d0/16.0d0*polter/ak2
c     TESTING
c        d=vis(j)*(0.25d0*deltag+alphadot)

        if (x.gt.0.0d0) then
           dp=vis(j)*3.0d0/16.0d0*polter/x**2
        else
           dp=0.0d0
        end if
cxxxc     TESTING
cxxx        if (a.gt.1.0d0/20.0d0) then
cxxx           d=vis(j)*(0.25d0*deltag
cxxx     2          +s2+polter/16+thetabdot/ak2+3.0d0/16.0d0
cxxx     3          *polterddot/ak2)+dvis(j)*(alpha+thetab/ak2+3.0d0/8.0d0
cxxx     4          *polterdot/ak2)+ddvis(j)*3.0d0/16.0d0*polter/ak2
cxxx           dp=expmmu(j)*s1
cxxx        else
cxxx           d=expmmu(j)*s1+vis(j)*(0.25d0*deltag
cxxx     2          +s2+polter/16+thetabdot/ak2+3.0d0/16.0d0
cxxx     3          *polterddot/ak2)+dvis(j)*(alpha+thetab/ak2+3.0d0/8.0d0
cxxx     4          *polterdot/ak2)+ddvis(j)*3.0d0/16.0d0*polter/ak2
cxxx           dp=0.0d0
cxxx        end if

c TESTING Newtonian gauge variables, Bardeen's curvature
c       psi=alphadot+adotoa*alpha
c        phi=eta-adotoa*alpha
c       phidot=etadot-adotoa*alphadot-adotoadot*alpha
c ignore massive neutrinos
c       grho=grhom*(omegac+omegab)/a+(grhog+grhor
c     1         *annur)/a2+grhom*omegavdyn*a2
c       gpres=((grhog+grhor*annur)/3.0d0)/a2
c     2     +weos*grhom*omegavdyn*a2
c       wprho=1+gpres/grho
c       zeta=-phi-2*(phidot/adotoa+psi)/3/wprho
c       deltamn=deltam-3*adotoa*alpha
c       deltabn=deltab-3*adotoa*alpha
c       write(16,'(8E13.4)')ak,a,(0.25d0*deltag+alphadot),
c     1  phi,psi,zeta,deltamn,deltabn

c lensing visibility function; approximate epoch of recombination
        phi=eta-adotoa*alpha
        chi=tau0-tau
        chir=tau0-taurmax
        if (chi.lt.chir) then
c lensing convergence, expmmu is an approximation
c one should integrate visibility function over r(chi)/r(tau)
c but the error is harmless
            glens=(chir-chi)*chi/chir
            dk=glens*ak2*phi*expmmu(j)
c alternative form; add l(l+1) in the final Cl
c it is numerically inaccurate at high l, equivalent to above at low l
c            glens=(chir-chi)/chi/chir
c            dk=glens*phi*expmmu(j)
c using Poisson's eq. below gives identical results
c            deltac=y(4)
c            dk=glens*grhom*(omegac*deltac
c     +omegab*deltab)/2.0/a*expmmu(j)
        else
           dk=0.0d0
        end if

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fderivst(n,x,y,yprime)
c  Evaluate the time derivatives of the perturbations.
c
        implicit double precision (a-h,o-z)
        dimension y(n),yprime(n)
c
        include 'cmbfast.inc'

        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $       itflag,lmaxt
c
c        common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
c        common /reionization2/ j2ri1,nri,nri0
c
        parameter (nqmax0=15)
        parameter(lmx0=30)
        parameter (ep0=1.0d-2)
c
        common /epsilon/ epsw
c
c  Internal variables.
        dimension denl(lmx0),dlfdlq(nqmax0)
        common /store/ denl,dlfdlq
c

c ep is used to stop the tight coupling approximation.
        if (ak.gt.0.06d0*epsw) then
           ep=ep0
        else
           ep=1.17d0*ep0
        end if

        tau=x
        a=y(1)
c
        a2=a*a
        call thermo(tau,cs2,opac,dlxedla)
c
c Tight Coupling parameters
        tcp=0.0d0
        tcp1=ak/opac
        tcp2=1.0d0/(opac*tau)
        if ((tcp1.gt.ep).or.(tcp2.gt.ep)) then
           tcp=1.0d0
        end if
c
c  Compute expansion rate.
        if (amnu.eq.0) then
           rhonu=1.0d0
           pnu=1.0d0/3.0d0
        else
           call nu1(a,rhonu,pnu)
        end if
        omegavdyn = omegav * dynrho(a)
c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
        grho=grhom*(omegac+omegab)/a+(grhog+grhor
     1       *annur+grhonr*annunr*rhonu)/a2
     2      +grhom*omegavdyn*a2

#ifdef DIM
c energy density for  5 dimensions
        grho=(sqrt(grho)+sqrt(omegav*grhom*a2))**2

#endif
        adotoa=sqrt(grho/3.0d0)
        yprime(1)=adotoa*a

        ind1=4
        ind2=ind1+lmaxt+1
        ind3=ind2+lmaxt+1
        shearg=y(ind1)/15+y(ind1+2)/21+y(ind1+4)/35
        shearr=y(ind3)/15+y(ind3+2)/21+y(ind3+4)/35
c only do relativistic components
        pi=(grhog*shearg+annur*grhor*shearr)/a2/3

c Tensors
        ht=y(2)
        htpr=y(3)
        yprime(2)=htpr
        htdpr=-2*adotoa*htpr-ak2*ht+24*pi
        yprime(3)=htdpr
c
c Photon perturbations
        psie=y(ind1)/10.0d0+y(ind1+2)/7.0d0+3d0*y(ind1+4)/70.0d0
     &       -3.0d0*y(ind2)/5.0d0+
     &      6.0d0*y(ind2+2)/7.0d0-3.0d0*y(ind2+4)/70.0d0
        if (tcp.eq.1) then
c no tight coupling approx
           yprime(ind1)=-ak*y(ind1+1)-opac*y(ind1)+opac*psie-htpr
           yprime(ind2)=-ak*y(ind2+1)-opac*y(ind2)-opac*psie
c l=1...lmaxt
           do 70 l=1,lmaxt-1
              yprime(ind1+l)=ak*denl(l)*(l*y(ind1-1+l)-(l+1)
     &   *y(ind1+1+l))-opac*y(ind1+l)
              yprime(ind2+l)=ak*denl(l)*(l*y(ind2-1+l)-(l+1)
     &   *y(ind2+1+l))-opac*y(ind2+l)
 70        continue
c
c Truncate moment expansion
           yprime(ind1+lmaxt)=ak*y(ind1-1+lmaxt)-(lmaxt+1)/tau
     2     *y(ind1+lmaxt)-opac*y(ind1+lmaxt)
           yprime(ind2+lmaxt)=ak*y(ind2-1+lmaxt)-(lmaxt+1)/tau
     2     *y(ind2+lmaxt)-opac*y(ind2+lmaxt)
        else
           deltat0=-4.0d0*htpr/opac/3.0d0
           deltap0=-deltat0/4.0d0
           y(ind1)=deltat0
           y(ind2)=deltap0
           do 80 l=0,lmaxt
              yprime(ind1+l)=0.0d0
              yprime(ind2+l)=0.0d0
 80        continue
        end if
c massless neutrinos (ignore massive)
           yprime(ind3)=-ak*y(ind3+1)-htpr
c l=1...lmaxt
           do 90 l=1,lmaxnr-1
              yprime(ind3+l)=ak*denl(l)*(l*y(ind3-1+l)-(l+1)
     &        *y(ind3+1+l))
 90        continue
c
c Truncate moment expansion
           yprime(ind3+lmaxnr)=ak*y(ind3-1+lmaxnr)-(lmaxnr+1)/tau
     2     *y(ind3+lmaxnr)
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine foutputt(n,y,ypr,j,tau0,tau,dt,dte,dtb)

        implicit double precision (a-h,o-z)
        integer j,n

        include 'cmbfast.inc'

        common /cosmoparm/ ak,ak2,amnu,lmaxg,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     $       itflag,lmaxt
c
        dimension y(n),ypr(n)
        dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
        dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
        common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
        save /visib/
c
c Tensors
c
        x=ak*(tau0-tau)
        x2=x*x
        ind1=4
        ind2=ind1+lmaxt+1
        htpr=y(3)
        htdpr=ypr(3)
        psie=y(ind1)/10.0d0+y(ind1+2)/7.0d0+3.0d0
     & *y(ind1+4)/70.0d0-3.d0*y(ind2)/5.0d0+6.0d0*y(ind2+2)
     & /7.0d0-3.0d0*y(ind2+4)/70.0d0
        psiedot=ypr(ind1)/10.0d0+ypr(ind1+2)/7.0d0+3.0d0
     & *ypr(ind1+4)/70.0d0-3.d0*ypr(ind2)/5.0d0+6.0d0*ypr(ind2+2)
     & /7.0d0-3.0d0*ypr(ind2+4)/70.0d0
        psieddot=-0.3d0*(opac(j)*psiedot+dopac(j)*psie)
     & -0.1d0*htdpr-ak*(3.0d0*ypr(ind1+1)/70.0d0+ypr(ind1+3)/15.0d0
     & +ypr(ind1+5)/42.0d0-33.0d0*ypr(ind2+1)/35.0d0
     & +8.0d0*ypr(ind2+3)/15.0d0-ypr(ind2+5)/42.0d0)
c

        if (x.gt.0.0d0) then
           dt=(-expmmu(j)*htpr+vis(j)*psie)/x2

           dte=vis(j)*(psie-psieddot/ak2-6.0d0*psie/x2
     &     -4.0d0*psiedot/ak/x)-dvis(j)*(4.0d0*psie/x/ak
     &     +2.0d0*psiedot/ak2)-ddvis(j)*psie/ak2
           dtb=2.0d0*(vis(j)*(2.0d0*psie/x+psiedot/ak)
     &     +dvis(j)*psie/ak)
        else
           dt=0.0d0
           dte=0.0d0
           dtb=0.0d0
        end if
c

        dte=-dte
        dtb=-dtb
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine finithermo(taumin,taumax,tau0,taurend,
     2     dlntau0,n1,nstep)
c     Compute and save unperturbed baryon temperature and ionization fraction
c     as a function of time.  With nthermo=10000, xe(tau) has a relative
c     accuracy (numerical integration precision) better than 1.e-5.
c
      implicit double precision (a-h,o-z)
c
      include 'cmbfast.inc'

      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr
      common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $     ,grhonr
c
      parameter (barssc0=9.1820d-14)
c      parameter (nstep0=2400)
c
      double precision atemp(nthermo),datemp(nthermo)
      dimension tb(nthermo),cs2(nthermo),xe(nthermo)
      dimension dcs2(nthermo)
      dimension dotmu(nthermo),sdotmu(nthermo),emmu(nthermo)
      dimension demmu(nthermo),ddotmu(nthermo)
      dimension dddotmu(nthermo),ddddotmu(nthermo)
      dimension vis(nstep0),dvis(nstep0),ddvis(nstep0)
      dimension opac(nstep0),dopac(nstep0),expmmu(nstep0)
      double precision atau0(nstep0),dtau1(nstep0),dtau2(nstep0)
      common /thermod/ tauminn,dlntau,dotmu,ddotmu,cs2,dcs2

      dimension dlxedla(nthermo),ddlxedla(nthermo)
      common /thermod2/ dlxedla,ddlxedla
      save /thermod2/

      common /visib/ vis,dvis,ddvis,expmmu,opac,dopac
      common /timesteps/ atau0,dtau1,dtau2
      common /reionization/zri,taurist,zristp,tauristp,rif,optdlss
      common /reionization2/ j2ri1,nri,nri0
      common /maxreion/ taurmax,armax,jrmax
      save /thermod/
      save /visib/
      integer rcflag
      common /recfl/rcflag

cxxxc TESTING xe
cxxx      character*500 filename
cxxx      parameter (nmax=nthermo)
cxxx      double precision xvec(nmax),yvec(nmax),dyvec(nmax)
cxxx


      ncount=0
      thomc0=5.0577d-8*tcmb**4
      akthom=(2.3048d-9)*(1-yhe)*omegab*h0*h0
      tauminn=0.05d0*taumin
      dlntau=log(tau0/tauminn)/(nthermo-1)
c
c     Initial conditions: assume radiation-dominated universe.
      tau01=tauminn
      adot0=adotrad
      astart=adotrad*tauminn
      a0=astart
      a02=a0*a0
c     Assume that any entropy generation occurs before tauminn.
c     This gives wrong temperature before pair annihilation, but
c     the error is harmless.
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
c     Integrate Friedmann equation using inverse trapezoidal rule.
         a=a0+adot0*dtau
         a2=a*a
         call nu1(a,rhonu,pnu)

         omegavdyn = omegav * dynrho(a)

         grho=grhom*(omegac+omegab)/a+(grhog+grhor
     1        *annur+grhonr*annunr*rhonu)/a2
     2        +grhom*omegavdyn*a2

#ifdef DIM
         grho=(sqrt(grhom*(omegac+omegab)/a+(grhog+grhor
     1          *(annur+annunr*rhonu))/a2
     2          +grhom*omegav*a2)
     3          +sqrt(omegav*grhom)*a)**2
#endif

         adot=sqrt(grho/3.0d0)*a
         a=a0+2.0d0*dtau/(1.0d0/adot0+1.0d0/adot)
         atemp(i)=a
c     Baryon temperature evolution: adiabatic except for Thomson cooling.
c     Use  quadrature solution.
         tg0=tcmb/a0
         ahalf=0.5d0*(a0+a)
         adothalf=0.5d0*(adot0+adot)
c     fe=number of free electrons divided by total number of free baryon
c     particles (e+p+H+He).  Evaluate at timstep i-1 for convenience; if
c     more accuracy is required (unlikely) then this can be iterated with
c     the solution of the ionization equation.
         fe=(1.d0-yhe)*xe(i-1)/(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i-1))
         thomc=thomc0*fe/adothalf/ahalf**3
         etc=exp(-thomc*(a-a0))
         a2t=a0*a0*(tb(i-1)-tg0)*etc-tcmb/thomc*(1.d0-etc)
         tb(i)=tcmb/a+a2t/(a*a)
c     Integrate ionization equation.
         tbhalf=0.5d0*(tb(i-1)+tb(i))

c     If there is re-ionization, smoothly increase xe to the
c     requested value.
         if ((zri.ne.0.0).and.(tau.gt.(9.0d0*taurist/10.0d0))) then
            if(ncount.eq.0) then
               ncount=i-1
            end if
c     SMOOTH REIONIZATION
            xod=150.0d0*(tau-taurist)/taurist
            if (xod.gt.100) then
               tgh=1
            else
               tgh=(exp(xod)-exp(-xod))/(exp(xod)+exp(-xod))
            end if
            xe(i)=(rif-xe(ncount))*(tgh+1.0d0)/2.0d0+xe(ncount)
            dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)
         else
            if (rcflag.eq.0) then
               call ionize(tbhalf,ahalf,adothalf,dtau,xe0)
               call ionhe(tb(i),a,xe0,x1,x2)
               xe(i)=xe0+0.25d0*yhe/(1.0d0-yhe)*(x1+2*x2)
               dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)
            else
               call recint(a,xe(i))

cxxxc     TESTING
cxxx               call splintrec(xvec,yvec,dyvec,nvec,log(a),xetemp)
cxxx               xetemp2=exp(xetemp)
cxxx               write(66,*)a,xe(i),xetemp2,
cxxx     $              0.5*(xetemp2-xe(i))/(xetemp2+xe(i))
cxxx               xe(i)=xetemp2

               dlxedla(i)=(xe(i)-xe(i-1))/dtau*a/adot/xe(i)

            endif
         end if

cxxxc     TESTING
cxxx         photbar=grhog/(grhom*omegab*a)
cxxx         pb43=4.0d0/3.0d0*photbar
cxxx         write(69,'(99E15.5)')tau,a,xe(i),dlxedla(i),
cxxx     $        2.0d0*pb43/(1.0d0+pb43)


c     Baryon sound speed squared (over c**2).
         dtbdla=-2.0d0*tb(i)-thomc*adothalf/adot*(a*tb(i)-tcmb)
         barssc=barssc0*(1.d0-0.75d0*yhe+(1.d0-yhe)*xe(i))
c         cs2(i)=barssc*tb(i)*(1-dtbdla/tb(i)/3.0d0)
         cs2(i)=5.0/3.0*barssc*tb(i)
c     Calculation of the visibility function
         dotmu(i)=xe(i)*akthom/a2
c
 15      a0=a
         tau01=tau
         adot0=adot
 10   continue

      if ((xe(nthermo).lt.rif).and.(zri.ne.0.0d0)) then
         write(*,*)'Warning: We use a smooth function to '
         write(*,*)'approach your specified reionization'
         write(*,*)'fraction. The redshift that is deduced from'
         write(*,*)'youre input paprameters is so low that our'
         write(*,*)'smooth function does not reach the required'
         write(*,*)'value. You should go in to subroutine finithermo'
         write(*,*)'and play with the shape of this smooth function.'
         write(*,*)'Search for SMOOTH REIONIZATION for'
         write(*,*)' the place where the function is set.'
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
      emmu(nthermo) = 0

      iv=0
      vfi=0.0d0
c     Getting the starting and finishing times for decoupling.
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
            taurst=9.0d0/10.0d0*tauminn*exp(dble(j1-1)*dlntau)
            iv=1
         end if
         if ((iv.eq.1).and.(vfi.gt.0.999)) then
            taurend1=(tauminn*exp(dble(j1-1)*dlntau))
            taurend1=max(taurend1,taurend1*sqrt(2500.0d0/
     &           (omegac+omegab)/h0**2))
            iv=2
         end if
 30   continue
      if (iv.ne.2) then
         taurend1=2.5d0*(tauminn*exp(dble(ncount-1)*dlntau))
      end if
c     Calculating the timesteps during recombination.

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

      if (dtaurec.ne.0.0) then
         dtaurec=min(dtaurec,taurst/40.0d0/accstep)
      else
         dtaurec=taurst/40.0d0/accstep
      end if

      taurend2=dtaurec/dlntau0
      taurend=max(taurend1,taurend2)
      taurend=min(taurend,9.0d0*taurist/10.0d0)
c     In models where reionization starts very early so that
c     it cut into what we call recombination we will
c     make the timesteps there at least as small as the
c     ones we were using during recombination. If not
c     timesteps 1.5 times bigger suffice.
      fact=1.5d0*accstep
      if (taurend.eq.(9.0d0*taurist/10.0d0))fact=1.0d0*accstep
      n1=int((taurend-taurst)/dtaurec+1)
      dtaurec=(taurend-taurst)/(n1-1)
      n1=n1+1
c     Calculating the timesteps after recombination (logarithmic
c     outside re-ionization scattering surface).
      nstep=int(log(taumax/taurend)/dlntau0)
      dlntau0=log(taumax/taurend)/nstep
      nstep=nstep+n1


c     Adjusting if there is reionization.
c     There will be nri0 points to sample the quick rise in
c     the free electron density. After that, timesteps of length
c     dtauri until tauristp.
c
      nri0=50
      if (zri.ne.0) then
c     modified for lensing
c         taurist1=10.0d0/9.0d0*taurend
         taurist1=taurist
         j2ri1=int(log(9.0d0*taurist1/taurend/10.0d0)
     &        /dlntau0+n1)

         tautemp=min(21.0d0*taurist1/20.0d0,tauristp)
         j2ri2=int(log(tautemp/taurend)
     &        /dlntau0+n1)

         j2ri3=int(log(tauristp/taurend)/dlntau0+n1)
         dtri0=taurend*(exp(dlntau0*(j2ri2-n1))-
     2        exp(dlntau0*(j2ri1-n1)))
         dtri=taurend*(exp(dlntau0*(j2ri3-n1))-
     2        exp(dlntau0*(j2ri2-n1)))
         dtauri0=dtri0/dble(nri0)

         dtauri=dtaurec*fact
         dtauri=min(dtauri,dtri/10.0d0)

c     FIXING
         if (dtauri.gt.0.0d0) then
            nri=int(dtri/dtauri)+1
            dtauri=dtri/dble(nri)
         else
            nri=0
            dtauri=0.0d0
         end if
         nstep=nstep+nri0+nri+j2ri1-j2ri3
      else
         j2ri1=0
         j2ri2=0
      end if


c
      if (nstep.gt.nstep0) then
         write(*,*)
     2        'Sorry, the arrays were dimensioned for a maximum of'
         write(*,*)nstep0, 'timesteps.'
         write(*,*)'The model you requested needs'
         write(*,*)nstep,'Please make the arrays bigger by making'
         write(*,*)'nstep0 bigger where it appears'
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
c         write(68,'(3E15.5)')tau,tvis,atemp(i)
c      end do
cc
c     Saving tau values and the quantities needed to calculate
c     the derivatives of the visibility function appearing in the sources.

c     TIMING
c        write(*,*)'timesteps',n1,j2ri1,nri
c        write(*,*)'timesteps',nri0,nstep
c        write(*,*)dtauri0,dtauri
c        write(*,*)taurend*exp(dlntau0*dble(j2ri1-n1)),tauristp
c        write(*,*)taurend*exp(dlntau0*dble(j2ri2-n1))
c        write(*,*)taurend*exp(dlntau0*dble(j2ri3-n1))
c        write(50,*)'timesteps',n1,j2ri1,nri
c        write(50,*)'timesteps',nri0,nstep

      vismax=0.0
      do 40 j2=2,nstep
         if (j2.le.n1) then
            tau=taurst+dble(j2-2)*dtaurec
         else
            if ((zri.eq.0).or.(j2.le.j2ri1)) then
               tau=taurend*exp(dlntau0*dble(j2-n1))
            else
               if (j2.lt.(j2ri1+nri+nri0)) then
                  if (j2.le.(j2ri1+nri0)) then
                     tau=atau0(j2ri1)+dtauri0*dble(j2-j2ri1)
                  else
                     tau=atau0(j2ri1+nri0)+dtauri
     &                    *dble(j2-j2ri1-nri0)
                  end if
               else
                  tau=taurend*exp(dlntau0*dble(j2-j2ri1-nri0
     2                 -nri+j2ri3-n1))
               end if
            end if
         end if
         atau0(j2)=tau
c     Cubic-spline interpolation.
         d=log(tau/tauminn)/dlntau+1.0d0
         i=int(d)
         d=d-i
         if (i.lt.nthermo) then
            opac(j2)=dotmu(i)+d*(ddotmu(i)
     $           +d*(3.0d0*(dotmu(i+1)-dotmu(i))
     2           -2.0d0*ddotmu(i)-ddotmu(i+1)
     $           +d*(ddotmu(i)+ddotmu(i+1)
     3           +2.0d0*(dotmu(i)-dotmu(i+1)))))
            dopac(j2)=(ddotmu(i)+d*(dddotmu(i)+d*(3.0d0*(ddotmu(i+1)
     2           -ddotmu(i))-2.0d0*dddotmu(i)
     $           -dddotmu(i+1)+d*(dddotmu(i)
     3           +dddotmu(i+1)+2.0d0*(ddotmu(i)
     $           -ddotmu(i+1))))))/(tau
     4           *dlntau)
            ddopac=(dddotmu(i)+d*(ddddotmu(i)
     $           +d*(3.0d0*(dddotmu(i+1)
     2           -dddotmu(i))-2.0d0*ddddotmu(i)-ddddotmu(i+1)
     3           +d*(ddddotmu(i)+ddddotmu(i+1)
     $           +2.0d0*(dddotmu(i)
     4           -dddotmu(i+1)))))-(dlntau**2)*tau*dopac(j2))
     5           /(tau*dlntau)**2
            expmmu(j2)=emmu(i)+d*(demmu(i)
     $           +d*(3.0d0*(emmu(i+1)-emmu(i))
     2           -2.0d0*demmu(i)-demmu(i+1)+d*(demmu(i)+demmu(i+1)
     3           +2.0d0*(emmu(i)-emmu(i+1)))))

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
     $           +3*opac(j2)*dopac(j2)+ddopac)
         else
            opac(j2)=dotmu(nthermo)
            dopac(j2)=ddotmu(nthermo)
            ddopac=dddotmu(nthermo)
            expmmu(j2)=emmu(nthermo)
            vis(j2)=opac(j2)*expmmu(j2)
            dvis(j2)=expmmu(j2)*(opac(j2)**2+dopac(j2))
            ddvis(j2)=expmmu(j2)*(opac(j2)**3+3.0d0*opac(j2)
     2           *dopac(j2)+ddopac)
c
         end if
 40   continue
      atau0(1)=0.0d0
      atau0(nstep)=min(tau0,atau0(nstep))
c     saving the length of the timesteps for the time integration.
      do 50 j2=2,nstep-1
         dtau1(j2)=atau0(j2+1)-atau0(j2)
         dtau2(j2)=abs(atau0(j2+1)-atau0(j2-1))/2.0d0
 50   continue

c     TESTING. Check visibility function.
c      write(*,*)'Maximum of vis'
c      write(*,*)taurmax,armax
cxxx      do j2=1,nstep
cxxx         write(60,'(1I6,4E15.5)')j2,atau0(j2),vis(j2)
cxxx      end do

      dtau1(nstep)=dtau1(nstep-1)/2.0d0
      dtau2(nstep)=(atau0(nstep)-atau0(nstep-1))/2.0d0
      dtau1(1)=(atau0(2)-atau0(1))/2.0d0
      dtau2(1)=dtau2(2)/2.0d0


c      stop

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine finitial(y,tau)
c  Initial conditions.
        implicit double precision (a-h,o-z)
        integer initfl
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                     nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /initcase/ initfl
        common /qparm/ wdyn,ndyn
c
        parameter (lmax0=12,lmaxnr0=25,lmaxnu0=25,nqmax0=15)
        parameter (nvar0=9+2*(lmax0+1)+(lmaxnr0+1)
     $       +nqmax0*(lmaxnu0+1))
        dimension y(nvar0)
c
        a=tau*adotrad
        a2=a*a
        call nu1(a,rhonu,pnu)
        omegavdyn = omegav * dynrho(a)
        weos = wdyn_func(a)

c  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
        grho=grhom*(omegac+omegab)/a+(grhog+grhor
     $       *annur+grhonr*annunr*rhonu)/a2
     2      +grhom*omegavdyn*a2

#ifdef DIM
c energy density for  5 dimensions
        grho=(sqrt(grho)+sqrt(omegav*grhom*a2))**2

#endif

c       adotoa=sqrt(grho/3.0d0)
        gpres=((grhog+grhor*annur)/3.0d0+grhonr*annunr*pnu)/a2
     2       +weos*grhom*omegavdyn*a2
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
c normalize to zeta=1
#ifdef OLDNORM
        psi=-1.0d0
#else
        psi=-1.0d0/(1.5d0+2.0d0*fracnu/5.0d0)
#endif
        C=(15.0d0+4.0d0*fracnu)/20.0d0*psi
        akt2=(ak*tau)**2
        h=C*akt2*(1.0d0-0.2d0*yrad)
        eta=2.0d0*C-(5.0d0+4.0d0*fracnu)/6.0d0/(15.0d0+4.0d0*fracnu)*
     2              C*akt2*(1.0d0-yrad/3.0d0)
        f1=(23.0d0+4.0d0*fracnu)/(15.0d0+4.0d0*fracnu)
        deltac=-0.5d0*h
        deltag=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0)
        deltab=0.75d0*deltag
        deltar=-2.0d0/3.0d0*h*(1.0d0-akt2/36.0d0*f1)
        thetac=0.0d0
        thetag=-C/18.0d0*akt2*akt2/tau
        thetab=thetag
        thetar=f1*thetag
        shearr=4.0d0/15.0d0*ak2/s*psi*(1.0d0+7.0d0/36.0d0*yrad)
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
        dgrho=grhom*(omegac*deltac+omegab*deltab)/a
     2    +(grhog*deltag+grhor*annur*deltar+grhonr*annunr*drhonu)/a2

c no initial perturbations in dark energy

       if (ndyn.eq.1.or.ndyn.eq.2) then
        rc_phi = 0.0
        rc_psi = 0.0
        y(nvar0-1) = rc_phi
        y(nvar0)   = rc_psi
         call dyn_nrg(a,grho,gpres,rc_phi,rc_psi,
     2          dgrho_phi,dgprs_phi,dgtheta_phi)
         dgrho=dgrho+dgrho_phi

       endif

c  Add a seed if desired.
        if (initfl.eq.4) dgrho=dgrho+grhom/a
c       econ=(adotoa*ahdot/a-2.0d0*ak2*eta-dgrho)/grho
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine finitialt(y,tau)
c  Initial conditions.
        implicit double precision (a-h,o-z)
        parameter (lmaxnr0=25)
c
        include 'cmbfast.inc'

        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec
     $       ,grhonr
        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax)
     $     ,itflag,lmaxt
c
        parameter (lmaxt0=10)
        parameter (nvar0t=2*(lmaxt0+1)+2+1+(lmaxnr0+1))
c
        dimension y(nvar0t)
c
        a=tau*adotrad
c
        y(1)=a
c Tensor modes
c       ht=1.0d0
c       htpr=0.0d0
        y(2)=1.0d0
        y(3)=0.0d0
        ind1=4
        ind2=ind1+lmaxt+1
        ind3=ind2+lmaxt+1
        do l=0,lmaxt
           y(ind1+l)=0.0d0
           y(ind2+l)=0.0d0
        end do
        do l=0,lmaxnr0
           y(ind3+l)=0.0d0
        end do
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initjl(filename)

c     This subroutine reads the jl files from disk and
c     initializes other variables needed in CMBFAST.

      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      parameter (ketamax=3*l0max+126)
      parameter (d0hi=1.0d40,d0lo=1.0d40,xlimmin=10.0d0)

      double precision x
      real aj
      double precision ajl(ketamax,lmax),ajlpr(ketamax,lmax)
      double precision xx(ketamax), dxx(ketamax)
      integer l(lmax),l0
      integer l0file,lfile(lmax)
      character*500 filename

      common /lvalues1/ l,l0,lmo
      common /lvalues2/ akmax0
      save /lvalues1/
      save /lvalues2/

      common /jlgen/ ajl,ajlpr,xx,dxx,kmax
      save /jlgen/

      open(unit=11,file=filename
     &             ,status='unknown',form='unformatted')
      rewind 11
      read(11)lmofile
      read(11)kmaxfile

c     TESTING
c      write(*,*)kmaxfile,'kmaxfile'

c The code needs 300 more ls for the lensing calculation.

      if ((lmo.gt.lmofile).or.(akmax0.gt.kmaxfile)) then
         write(*,*)'You have entered a lmax and/or kmax'
         write(*,*)'inconsistent with those in the file'
         write(*,*)lmofile-300,kmaxfile
         write(*,*)'You will have to start again'
         stop
      end if

      kmaxfile=kmaxfile-25+151
      kmax=int(akmax0-25+151)

      if (kmaxfile.gt.ketamax) then
         write(*,*)'kmax in file =', kmaxfile,'is too large'
         write(*,*)'ketamax in routines is dimentioned to'
         write(*,*) ketamax
         write(*,*)'You will have to increase l0max'
         write(*,*)'in cmbfast.inc to increase ketamax=3*l0max+126'
         write(*,*)'it consistently everywhere'
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
c     reading  j_l
c     remember to create jl.dat with jlgen.f first using
c     correct lmax and ketamax

      do 40 i=1,kmaxfile
         if (i.le.151) then
            if (i.le.51) then
               xx(i)=dble((i-1))/10.0d0
            else
               xx(i)=dble((i-51))/5.0d0+5.0d0
            end if
         else
            xx(i)=dble((i-151))+25.0d0
         end if
 40   continue
      do 50 j=1,l0
         do 60 i=1,kmaxfile
            x=xx(i)
            xlim=dble(l(j))/20.0d0
            xlim=max(xlim,xlimmin)
            xlim=l(j)-xlim
c            if (x.gt.xlim) then
            if (x-xlim.gt.1E-10) then
               read(11)aj
               ajl(i,j)=dble(aj)
            else
               ajl(i,j)=0.0d0
            end if
 60      continue
 50   continue
      close(11)
      do 70 i=2,(kmaxfile-1)
         dxx(i)=(xx(i+1)-xx(i-1))/2.0d0
 70   continue
      dxx(1)=xx(2)/2.0d0
      dxx(kmaxfile)=(xx(kmaxfile)-xx(kmaxfile-1))/2.0d0
c
c     get the interpolation matrix for bessel functions

      do 80 j=1,l0
         call spline(xx,ajl(1,j),kmaxfile,d0lo,d0hi,ajlpr(1,j))
 80   continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine powersflat(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the scalar perturbations.

      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr


c     K split
      common /cutk/ aksplit, kcutflag
      if (kcutflag.eq.0) then
         win=1.0d0
      else
         win=2.0d0*(ak/aksplit)**4
         win=2.0d0*exp(-win)/(1.0d0+exp(-win))
         if (kcutflag.eq.-1) win=1.0d0-win
      end if

c     Normalize so tilt does not change power at
c     pivot point k=0.05/Mpc
      aknlog=log(ak/0.05d0)
      apower=exp((an(in)-1.0d0+.5*alphans(in)*aknlog+
     $     dalphansdlnk(in)*aknlog*aknlog/6.)*aknlog)


c     K split
      apower=apower*win

#ifdef OLDNORM
      apower=apower*(4.0d0*3.14159265d0)**2/2.0d0*100.0d0
#endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine powertflat(ak,in,apower)

c This subroutine computes the power spectra
c for mode ak of the tensor perturbations.

      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax)
     $     ,itflag,lmaxt
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr

#ifdef OLDNORM
      aknlog=log(ak/0.002d0)
#else
      aknlog=log(ak/0.05d0)
#endif
      apower=exp((ant(in)+.5*alphant(in)*aknlog)*aknlog)

#ifdef OLDNORM
      apower=apower*(4.0d0*3.14159265d0)**2*8.0d0*2.0d4
#endif
      return
      end

cxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxxxccccccc ADDED TO TEST RECOMBINATION ccccccccccccccccccccccccccccccccccc
cxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxxx
cxxx    subroutine makelogspline(filename,xvec,yvec,dyvec,nvec,i1,i2)
cxxx    implicit none
cxxx        integer nmax
cxxx        parameter (nmax=40000)
cxxx    double precision xvec(nmax),yvec(nmax),dyvec(nmax)
cxxx    integer nvec,i1,i2
cxxx    double precision temp(5)
cxxx    character*500 filename
cxxx    integer i
cxxx
cxxx    if((i1.gt.5).or.(i2.gt.5)) then
cxxx       write(*,*)'i1,i2 out of range'
cxxx       stop
cxxx    end if
cxxx
cxxx
cxxx    open(unit=10,file=filename,status='old',form='formatted')
cxxx
cxxx    i=0
cxxx 10 i=i+1
cxxx    if (i.gt.nmax) then
cxxx       write(*,*)'Arrays not big enough'
cxxx       stop
cxxx    end if
cxxx    read(10,*,end=20)temp(1),temp(2)
cxxx    xvec(i)=log(1.0d0/(temp(i1)+1.0d0))
cxxx    yvec(i)=log(temp(i2))
cxxx    goto 10
cxxx 20 continue
cxxx    close(10)
cxxx    nvec=i-1
cxxx
cxxx    call splinerec(xvec,yvec,nvec,1.0d40,1.0d40,dyvec)
cxxx
cxxx    return
cxxx    end
cxxxccccccccccccccccccccxsccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxxx      SUBROUTINE splintrec(xa,ya,y2a,n,x,y)
cxxx      INTEGER n
cxxx      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
cxxx      INTEGER k,khi,klo
cxxx      DOUBLE PRECISION a,b,h
cxxx      klo=1
cxxx      khi=n
cxxx1     if (khi-klo.gt.1) then
cxxx        k=(khi+klo)/2
cxxx        if(xa(k).gt.x)then
cxxx          khi=k
cxxx        else
cxxx          klo=k
cxxx        endif
cxxx      goto 1
cxxx      endif
cxxx      h=xa(khi)-xa(klo)
cxxx      if (h.eq.0.d0) pause 'bad xa input in splint'
cxxx      a=(xa(khi)-x)/h
cxxx      b=(x-xa(klo))/h
cxxx      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
cxxx     *2)/6.d0
cxxx      return
cxxxC  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
cxxx      END
cxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxxx
cxxx      SUBROUTINE splinerec(x,y,n,yp1,ypn,y2)
cxxx      INTEGER n,NMAX
cxxx      double precision yp1,ypn,x(n),y(n),y2(n)
cxxx      PARAMETER (NMAX=100010)
cxxx      INTEGER i,k
cxxx      double precision p,qn,sig,un,u(NMAX)
cxxx      if (yp1.gt..99d30) then
cxxx        y2(1)=0.d0
cxxx        u(1)=0.d0
cxxx      else
cxxx        y2(1)=-0.5d0
cxxx        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
cxxx      endif
cxxx      do 11 i=2,n-1
cxxx        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
cxxx        p=sig*y2(i-1)+2.d0
cxxx        y2(i)=(sig-1.d0)/p
cxxx        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
cxxx     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
cxxx     *u(i-1))/p
cxxx11    continue
cxxx      if (ypn.gt..99d30) then
cxxx        qn=0.d0
cxxx        un=0.d0
cxxx      else
cxxx        qn=0.5d0
cxxx        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
cxxx      endif
cxxx      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
cxxx      do 12 k=n-1,1,-1
cxxx        y2(k)=y2(k)*y2(k+1)+u(k)
cxxx12    continue
cxxx      return
cxxxC  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
cxxx      END
cxxx
cxxxCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
