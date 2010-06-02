C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (mzaldarriaga@cfa.harvard.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine output(clts,cltt,cles,clet,clbt,clcs,clct,
     &  itflag,lmx)

        implicit double precision(a-h,o-z)

        include 'cmbfast.inc'

        double precision clts(l0max,nnmax),cltt(l0max,nnmax)
        double precision cles(l0max,nnmax),clet(l0max,nnmax)
        double precision clbt(l0max,nnmax)
        double precision clcs(l0max,nnmax),clct(l0max,nnmax)


      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn

c Do not output last 300 ls, they are used internally to make a more
c accurate lensing calculation.

        lmxout=lmx-300

c Output
        do in=1,nn
           do il=2,lmxout
              if (itflag.ne.2) then
                 write(8,'(1I5,3E13.5)')il,clts(il,in)
     &           ,cles(il,in),clcs(il,in)
              end if
              if (itflag.ne.0) then
                 write(9,'(1I5,4E13.5)')il,cltt(il,in)
     &                      ,clet(il,in),clbt(il,in),clct(il,in)
              end if
           end do
        end do
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine COBEnormalize(clts,cltt,
     &  cles,clet,clbt,clcs,clct,clkk,cltk)

        implicit double precision(a-h,o-z)

        include 'cmbfast.inc'

        parameter (sigr0=8.0d0)
        parameter (fourpi=4.0d0*3.14159265d0)

        double precision ztf(ntfmax)
        double precision s8(ntfmax,nnmax),anorm(ntfmax,nnmax)

        double precision clts(l0max,nnmax),cltt(l0max,nnmax)
        double precision cles(l0max,nnmax),clet(l0max,nnmax)
        double precision clbt(l0max,nnmax)
        double precision clcs(l0max,nnmax),clct(l0max,nnmax)
        double precision clkk(l0max,nnmax),cltk(l0max,nnmax)

        integer l(lmax),l0
        double precision xl(lmax)
        common /lvalues1/ l,l0,lmo
        save /lvalues1/

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn

        common /tensor/ant(nnmax),rat(nnmax),alphant(nnmax),
     1       itflag,lmaxt

        common /transfer/akmaxt,ztf,nlnkt,ict,ntf

        common /transfernorm/ s8,anorm

        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /lensnorm/ d2normlens

        xlog10=log(10.0d0)
        h=h0/100.0d0
        xlnh=log(h)

        do j=1,l0
           xl(j)=dble(l(j))
        end do
c Curvature radius
        if (abs(omegak).gt.1.0d-3) then
           hc=2.998d5/h0
           curv=-omegak/hc/hc
           r=1.0d0/sqrt(abs(curv))
        endif


c Ensuring scalar to tensor ratio
        if (itflag.eq.1) then
           do in=1,nn
              do il=2,l(l0)
                 cltt(il,in)=cltt(il,in)*rat(in)
                 clet(il,in)=clet(il,in)*rat(in)
                 clbt(il,in)=clbt(il,in)*rat(in)
                 clct(il,in)=clct(il,in)*rat(in)
                 clkk(il,in)=clkk(il,in)*rat(in)
                 cltk(il,in)=cltk(il,in)*rat(in)
              end do
           end do
        end if
c COBE normalization
c fit the spectrum to a quadratic around C_10 with equal weights in logl

        do 700 in=1,nn

           c10=clts(10,in)+cltt(10,in)

           d1=(clts(int(xl(2)),in)+cltt(int(xl(2)),in))/c10-1.0d0
           d2=(clts(int(xl(3)),in)+cltt(int(xl(3)),in))/c10-1.0d0
           d3=(clts(int(xl(5)),in)+cltt(int(xl(5)),in))/c10-1.0d0
           d4=(clts(int(xl(7)),in)+cltt(int(xl(7)),in))/c10-1.0d0
           d5=(clts(int(xl(10)),in)+cltt(int(xl(10)),in))/c10-1.0d0
           d6=(clts(int(xl(11)),in)+cltt(int(xl(11)),in))/c10-1.0d0
           d7=(clts(int(xl(12)),in)+cltt(int(xl(12)),in))/c10-1.0d0

           x1=log(xl(2))/xlog10-1.0d0
           x2=log(xl(3))/xlog10-1.0d0
           x3=log(xl(5))/xlog10-1.0d0
           x4=log(xl(7))/xlog10-1.0d0
           x5=log(xl(10))/xlog10-1.0d0
           x6=log(xl(11))/xlog10-1.0d0
           x7=log(xl(12))/xlog10-1.0d0
           sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7
           s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7
           sx=x1**3+x2**3+x3**3+x4**3+x5**3+x6**3+x7**3
           sxy=x1**2*d1+x2**2*d2+x3**2*d3+x4**2*d4+
     2        x5**2*d5+x6**2*d6+x7**2*d7
           sxx=x1**4+x2**4+x3**4+x4**4+x5**4+x6**4+x7**4
           delt=s*sxx-sx*sx
           d1pr=(sxx*sy-sx*sxy)/delt
           d1ppr=2.0d0*(s*sxy-sx*sy)/delt

c Bunn and White fitting formula
           c10=(0.64575d0+0.02282d0*d1pr+0.01391d0*d1pr*d1pr
     2     -0.01819d0*d1ppr-0.00646d0*d1pr*d1ppr
     2     +0.00103d0*d1ppr*d1ppr)/c10
c logl
c           xlogl=-0.01669d0+1.19895d0*d1pr-0.83527d0*d1pr*d1pr
c     2           -0.43541d0*d1ppr-0.03421d0*d1pr*d1ppr
c     2           +0.01049d0*d1ppr*d1ppr
c          write(*,*)'COBE Likelihood relative to flat=',exp(xlogl)

c density power spectrum normalization;

           if (ict.eq.2) then
            do itf=1,ntf
#ifdef OLDNORM
             d2norm=c10*1.1d-9/clts(1,in)
#else
             d2norm=c10*1.1d-9/clts(1,in)/fourpi*2
#endif

             rewind(12+itf)
             dsig8=0.0d0
             dsig8o=0.0d0
             sig8=0.0d0
             alnko=0.0d0
             ict0=0

22           continue
             if (omegan.eq.0.0) then
              read(12+itf,*,end=23)akh,tfc,tfb,tfg,tfn
              tfm=(tfc*omegac+tfb*omegab)/(omegac+omegab)
             else
              read(12+itf,*,end=23)akh,tfc,tfb,tfg,tfn,tfnr
              tfm=(tfc*omegac+tfb*omegab+tfnr*omegan)/
     1        (omegac+omegab+omegan)
             endif
             ak=akh*h
             x=akh*sigr0
             win=3.0d0*(sin(x)-x*cos(x))/(x*x*x)
             alnk=log(ak)
             dlnk=alnk-alnko
             if (ict0.eq.0) then
              dlnk=2.0d0/(an(in)+3.0d0)
              ict0=1
             endif
             if (abs(omegak).lt.1.0d-3) then
                call powersflat(ak,in,apowers)
             else
                call powersopen(ak,in,apowers)
             endif
             dsig8=win*win*exp(alnk*4.0d0)*tfm*tfm*apowers
             sig8=sig8+(dsig8+dsig8o)*dlnk/2.0d0
             dsig8o=dsig8
             alnko=alnk
             goto 22
23           continue
             sig8=sqrt(sig8*d2norm)

             write(*,*)'sigma8=',sig8,' 4pik^3P(k) norm:',d2norm
             d2normlens=d2norm
c             write(12+itf,*)'sigma8=',sig8,' 4pik^3P(k) norm:',
c     1         d2norm

c       Change d2norm to be consistent with lensing,
c       where we use k/h
c            d2norm=d2norm*exp(xlnh*3.0d0)
             d2norm=d2norm*exp(xlnh*4.0d0)
             s8(itf,in)=sig8
             anorm(itf,in)=d2norm/fourpi
            enddo
           endif

c C_l normalization; output l(l+1)C_l/twopi
           c10=c10*2.2d-9/fourpi
           do il=2,l(l0)
              clts(il,in)=clts(il,in)*c10
              cles(il,in)=cles(il,in)*c10
              clcs(il,in)=clcs(il,in)*c10
              cltt(il,in)=cltt(il,in)*c10
              clet(il,in)=clet(il,in)*c10
              clbt(il,in)=clbt(il,in)*c10
              clct(il,in)=clct(il,in)*c10
              clkk(il,in)=clkk(il,in)*c10
              cltk(il,in)=cltk(il,in)*c10
           end do
 700       continue
           return
           end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine initlval(lmoin,akmax0in)

c This subroutines initializes l arrays.

        implicit double precision(a-h,o-z)

        include 'cmbfast.inc'

        integer l(lmax),l0

        common /lvalues1/ l,l0,lmo
        common /lvalues2/ akmax0
        save /lvalues1/
        save /lvalues2/

        common /lens/ lensflag


        akmax0=akmax0in
c 300 more ls will always be calculated, this is necessary
c for lensing.

        lmo=lmoin+300

        include 'lvalues.inc'

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine thermo(tau,cs2b,opac,dxe)
c  Compute unperturbed baryon temperature, sound speed squared,
c  and ionization fraction by interpolating pre-computed tables.
c
        implicit double precision (a-h,o-z)
c
        include 'cmbfast.inc'

        dimension cs2(nthermo),dcs2(nthermo)
        dimension dotmu(nthermo),ddotmu(nthermo)
        common /thermod/ tauminn,dlntau,dotmu,ddotmu,cs2,dcs2
        save /thermod/

c     TESTING
        dimension dlxedla(nthermo),ddlxedla(nthermo)
        common /thermod2/ dlxedla,ddlxedla
        save /thermod2/
c
        d=log(tau/tauminn)/dlntau+1.0d0
        i=int(d)
        d=d-i
        if (i.lt.1) then
c  Linear interpolation if out of bounds (should not occur).
          cs2b=cs2(1)+(d-1)*dcs2(1)
          opac=dotmu(1)+(d-1)*ddotmu(1)
          dxe=dlxedla(1)+(d-1)*ddlxedla(1)
        else if (i.ge.nthermo) then
          cs2b=cs2(nthermo)+(d-nthermo)*dcs2(nthermo)
          opac=dotmu(nthermo)+(d-nthermo)*ddotmu(nthermo)
          dxe=dlxedla(nthermo)+(d-nthermo)*ddlxedla(nthermo)
        else
c  Cubic spline interpolation.
          cs2b=cs2(i)+d*(dcs2(i)+d*(3.0d0*(cs2(i+1)-cs2(i))
     2         -2.0d0*dcs2(i)-dcs2(i+1)+d*(dcs2(i)+dcs2(i+1)
     3         +2.0d0*(cs2(i)-cs2(i+1)))))
          opac=dotmu(i)+d*(ddotmu(i)+d*(3.0d0*(dotmu(i+1)-dotmu(i))
     2         -2.0d0*ddotmu(i)-ddotmu(i+1)+d*(ddotmu(i)+ddotmu(i+1)
     3         +2.0d0*(dotmu(i)-dotmu(i+1)))))
          dxe=dlxedla(i)+d*(ddlxedla(i)+
     1         d*(3.0d0*(dlxedla(i+1)-dlxedla(i))
     2         -2.0d0*ddlxedla(i)-ddlxedla(i+1)+
     2         d*(ddlxedla(i)+ddlxedla(i+1)
     3         +2.0d0*(dlxedla(i)-dlxedla(i+1)))))

       end if
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine ionize(tempb,a,adot,dtau,xe)
c  Integrate the ionization fraction xe for hydrogen semi-implicitly
c  from tau to tau+dtau, treating tempb, a, and adot as constants
c  during this interval.
c
        implicit double precision (a-h,o-z)
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
c  Ionization temperature and coefficient.
        parameter (tion=1.5789d5,beta0=43.082d0)
c  Two-photon decay rate (in 1/Mpc).
        parameter (dec2g=8.468d14)
c  Switch for implicit (switch=1.0) or semi-implicit (switch=0.5) scheme.
        parameter (switch=0.5d0)
c
c  Recombination coefficient (in sqrt(K)/Mpc).
        alpha0=2.3866d-6*(1-yhe)*omegab*h0*h0
c  Coefficient for correction of radiative decay (dimensionless).
        crec=8.0138d-26*(1-yhe)*omegab*h0*h0
c  Recombination and ionization rates.
        phi2=0.448d0*log(tion/tempb)
        phi2=max(phi2,0.0d0)
        alpha=alpha0/sqrt(tempb)*phi2/(a*a*a)
        beta=tempb*phi2*exp(beta0-tion/tempb)
c  Peebles' correction factor.
        if (tempb.lt.200.0d0) then
          cpeebles=1.0d0
        else
          cp1=crec*dec2g*(1.0d0-xe)/(a*adot)
          cp2=crec*tempb*phi2*exp(beta0-0.25d0*tion/tempb)*
     2        (1.0d0-xe)/(a*adot)
          cpeebles=(1.0d0+cp1)/(1.0d0+cp1+cp2)
        end if
c  Integrate dxe=bb*(1-xe)-aa*xe*xe by averaging rhs at current tau
c  (fraction 1-switch) and future tau (fraction switch).
        aa=a*dtau*alpha*cpeebles
        bb=a*dtau*beta*cpeebles
        b1=1.0d0+switch*bb
        bbxe=bb+xe-(1.0d0-switch)*(bb*xe+aa*xe*xe)
        rat=switch*aa*bbxe/(b1*b1)
c  Prevent roundoff error.
        if (rat.lt.1.0d-6) then
          xe=bbxe/b1*(1.0d0-rat)
        else
          xe=b1/(2.0d0*switch*aa)*(sqrt(4.0d0*rat+1.0d0)-1.0d0)
        end if
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ionhe(tempb,a,x0,x1,x2)
c  Compute the helium ionization fractions using the Saha equation.
c  x0 is the hydrogen ionization fraction n(H+)/n(H) (input),
c  x1 is the helium first ionization fraction n(He+)/n(He)
c    (input and output), and
c  x2 is the helium second ionization fraction n(He++)/n(He)
c    (input and output).
c
        implicit double precision (a-h,o-z)
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
c  Ionization temperatures.
        parameter (tion1=2.855d5,tion2=6.313d5)
c
c  Constant for electron partition function per baryon.
        b0=2.150d24/((1.0d0-yhe)*omegab*h0*h0)
c  Electron partition function per baryon.
        b=b0*a*a*a*tempb*sqrt(tempb)
c  Dimensionless right-hand sides in Saha equations.
        if (abs(tion1/tempb).lt.100) then
           r1=4.0d0*b*exp(-tion1/tempb)
        else
           r1=0.0d0
        end if
        if (abs(tion2/tempb).lt.150) then
           r2=b*exp(-tion2/tempb)
        else
           r2=0.0d0
        end if
c
c  Solve coupled equations iteratively.
        c=0.25d0*yhe/(1.0d0-yhe)
        err=1.0d0
        niter=0
10      niter=niter+1
        if (err.lt.1.0d-12) return
          xe=x0+c*(x1+2*x2)
          x2new=r1*r2/(r1*r2+xe*r1+xe*xe)
          x1=xe*r1/(r1*r2+xe*r1+xe*xe)
          err=abs(x2new-x2)
          x2=x2new
          go to 10
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nu1(a,rhonu,pnu)
c  Compute massive neutrino density and pressure in units of the mean
c  density of one flavor of massless neutrinos.  Use cubic splines to
c  interpolate from a table.
c
        implicit double precision (a-h,o-z)
c
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                            nqmax,iq0,iq1,iq2
        parameter (nrhopn=10000)
        dimension r1(nrhopn),p1(nrhopn)
        dimension dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn)
        common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
        save /nu1d/

        if (amnu.eq.0.0d0) then
          rhonu=1.0d0
          pnu=1.0d0/3.0d0
          return
        end if

        d=log(a/amin)/dlna+1.0d0
        i=int(d)
        d=d-i
        if (i.lt.1) then
c  Use linear interpolation, bounded by results for massless neutrinos.
          rhonu=r1(1)+(d-1)*dr1(1)
          pnu=p1(1)+(d-1)*dp1(1)
          rhonu=min(exp(rhonu),1.0d0)
          pnu=min(exp(pnu),0.3333333333d0)
        else if (i.ge.nrhopn) then
c  This should not happen, unless the user evolves to z<0!
          rhonu=r1(nrhopn)+(d+i-nrhopn)*dr1(nrhopn)
          pnu=p1(nrhopn)+(d+i-nrhopn)*dp1(nrhopn)
          rhonu=exp(rhonu)
          pnu=exp(pnu)
        else
c  Cubic spline interpolation.
          rhonu=r1(i)+d*(dr1(i)+d*(3.0d0*(r1(i+1)-r1(i))
     1          -2.0d0*dr1(i)
     2         -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2.0d0*(r1(i)-r1(i+1)))))
          pnu=p1(i)+d*(dp1(i)+d*(3.0d0*(p1(i+1)-p1(i))-2.0d0*dp1(i)
     2         -dp1(i+1)+d*(dp1(i)+dp1(i+1)+2.0d0*(p1(i)-p1(i+1)))))
          rhonu=exp(rhonu)
          pnu=exp(pnu)
        end if

c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine initnu1(amnu)
c  Initialize interpolation tables for massive neutrinos.
c  Use cubic splines interpolation of log rhonu and pnu vs. log a.
        implicit double precision (a-h,o-z)
c
        parameter (nrhopn=10000,nqmax0=15)
        dimension r1(nrhopn),p1(nrhopn)
        dimension dr1(nrhopn),dp1(nrhopn)
        dimension ddr1(nrhopn)
        dimension qdn(nqmax0)
        common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
        save /nu1d/
        common /fermi/ qdn
        save /fermi/
c
        amin=1.0d-9
        dlna=-log(amin)/(nrhopn-1)
c
c  Check whether correct interpolation table already exists on disk.
        open(12,err=10,file='nu1.dat',status='old',form='unformatted')
        read(12,end=10) amnu1
        err=abs(amnu1/amnu-1.0d0)
        if (err.gt.1.0e-8) then
          go to 10
        else
          read(12,end=10) (r1(j),p1(j),j=1,nrhopn)
          close(12)
c  Yes, interpolation table was okay.
          go to 20
        end if

c  No, correct interpolation tables aren't on disk: compute and save.
10      close(12)
        do i=1,nrhopn
          a=amin*exp((i-1)*dlna)
          call ninu1(a,rhonu,pnu,amnu)
          r1(i)=log(rhonu)
          p1(i)=log(pnu)
        end do

20      call splini
        call splder(r1,dr1,nrhopn)
        call splder(p1,dp1,nrhopn)
        call splder(dr1,ddr1,nrhopn)

        open(12,file='nu1.dat',status='unknown',form='unformatted')
        rewind 12
        write(12) amnu
        write(12) (r1(j),p1(j),j=1,nrhopn)
        close(12)
c Calculate qdn
        dq=1.0d0
          do 30 i=1,nqmax0
          q=i-0.5d0
c         aq=a*amnu/q
c         v=1.0d0/sqrt(1.0d0+aq*aq)
          qdn(i)=dq*q*q*q/(exp(q)+1.0d0)
30      continue

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ninu1(a,rhonu,pnu,amnu)
c  Compute the density and pressure of one flavor of massive neutrinos,
c  in units of the mean density of one flavor of massless neutrinos.
c  Relative accuracy is better than 1.e-10 for all a.
c
        implicit double precision (a-h,o-z)
c
c       common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
c     &                              nqmax,iq0,iq1,iq2
c
        parameter (qmax=30.0d0,nq=1000,nq1=nq+1)
        dimension dum1(nq1),dum2(nq1)
c  const=7*pi**4/120.
        parameter (const=5.68219698d0)
c
        if (amnu.eq.0.0d0) then
          rhonu=1.0d0
          pnu=1.0d0/3.0d0
          return
        end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
c  Integrate up to qmax and then use asymptotic expansion for remainder.
        dq=qmax/nq
        dum1(1)=0.0d0
        dum2(1)=0.0d0
          do 10 i=1,nq
          q=i*dq
          aq=a*amnu/q
          v=1.0d0/sqrt(1.0d0+aq*aq)
          qdn=dq*q*q*q/(exp(q)+1.0d0)
          dum1(i+1)=qdn/v
          dum2(i+1)=qdn*v
10      continue
        call splint(dum1,rhonu,nq1)
        call splint(dum2,pnu,nq1)
c  Apply asymptotic corrrection for q>qmax and normalize by relativistic
c  energy density.
        rhonu=(rhonu+dum1(nq1)/dq)/const
        pnu=(pnu+dum2(nq1)/dq)/const/3.0d0
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nu2(a,drhonu,fnu,dpnu,shearnu,psi0,psi1,psi2)
c  Compute the perturbations of density, energy flux, pressure, and
c  shear stress of one flavor of massive neutrinos, in units of the mean
c  density of one flavor of massless neutrinos, by integrating over
c  momentum.
c
        implicit double precision (a-h,o-z)
c
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                        nqmax,iq0,iq1,iq2
        common /fermi/ qdn
        save /fermi/
c
        parameter (nqmax0=15,qmax=nqmax0-0.5d0)
        dimension psi0(nqmax0),psi1(nqmax0),psi2(nqmax0)
c  const=7*pi**4/120.
        parameter (const=5.68219698d0)
c
        dimension g0(4),g1(nqmax0+1),g2(nqmax0+1)
        dimension g3(nqmax0+1),g4(nqmax0+1)
        dimension qdn(nqmax0)
c
        if (nqmax.eq.0) then
          drhonu=0.0d0
          fnu=0.0d0
          dpnu=0.0d0
          shearnu=0.0d0
          return
        end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0.0d0
        g2(1)=0.0d0
        g3(1)=0.0d0
        g4(1)=0.0d0
        do 40 iq=2,(nqmax0+1)
            q=iq-1.5d0
            aq=a*amnu/q
            v=1.0d0/sqrt(1.0d0+aq*aq)
            g1(iq)=qdn(iq-1)*psi0(iq-1)/v
            g2(iq)=qdn(iq-1)*psi0(iq-1)*v
            g3(iq)=qdn(iq-1)*psi1(iq-1)
            g4(iq)=qdn(iq-1)*psi2(iq-1)*v
40      continue
        call splint(g1,g0(1),nqmax0+1)
        call splint(g2,g0(2),nqmax0+1)
        call splint(g3,g0(3),nqmax0+1)
        call splint(g4,g0(4),nqmax0+1)
        gf1=g1(nqmax0+1)
        gf2=g2(nqmax0+1)
        gf3=g3(nqmax0+1)
        gf4=g4(nqmax0+1)
        drhonu=(g0(1)+gf1*2.0d0/qmax)/const
        fnu=(g0(3)+gf3*2.0d0/qmax)/const
        dpnu=(g0(2)+gf2*2.0d0/qmax)/const/3.0d0
        shearnu=(g0(4)+gf4*2.0d0/qmax)/const*2.0d0/3.0d0
c
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine nuder(a,adotoa,rhonu,rhonudot,shearnudot,
     2                  psi2,psi2dot)
c  Compute the time derivative of the mean density in massive neutrinos
c  and the shear perturbation.
c
        implicit double precision (a-h,o-z)
c
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,
     &                    lmaxnu,nqmax,iq0,iq1,iq2
        common /fermi/ qdn
c
        parameter (nrhopn=10000)
        dimension r1(nrhopn),p1(nrhopn)
        dimension dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn)
        common /nu1d/ amin,dlna,r1,p1,dr1,dp1,ddr1
c
        parameter (nqmax0=15,qmax=nqmax0-0.5d0)
        dimension psi2(nqmax0),psi2dot(nqmax0)
c
        parameter (const=5.68219698d0)
c
        dimension g1(nqmax0+1)
        dimension qdn(nqmax0)
c
        save /nu1d/
        save /fermi/
c
        if (nqmax.eq.0) then
          rhonudot=0.0d0
          shearnudot=0.0d0
          return
        end if
c
c  q is the comoving momentum in units of k_B*T_nu0/c.
        g1(1)=0.0d0
        do 40 iq=2,(nqmax0+1)
            q=iq-1.5d0
            aq=a*amnu/q
            aqdot=aq*adotoa
            v=1.0d0/sqrt(1.0d0+aq*aq)
            vdot=-aq*aqdot/(1.0d0+aq*aq)**1.5d0
            g1(iq)=qdn(iq-1)*(psi2dot(iq-1)*v+psi2(iq-1)*vdot)
40      continue
        call splint(g1,g0,nqmax0+1)
        gf1=g1(nqmax0+1)
        shearnudot=(g0+gf1*2.0d0/qmax)/const*2.0d0/3.0d0
c
c
        d=log(a/amin)/dlna+1.0d0
        i=int(d)
        d=d-i
        if (i.lt.1) then
c  Use linear interpolation
          rhonudot=dr1(1)+(d-1)*ddr1(1)
        else if (i.gt.nrhopn) then
c  This should not happen, unless the user evolves to z<0!
          rhonudot=dr1(nrhopn)+(d+i-nrhopn)*ddr1(nrhopn)
        else

c  Cubic spline interpolation for rhonudot.
          rhonudot=dr1(i)+d*(ddr1(i)+d*(3.0d0*(dr1(i+1)-dr1(i))
     2         -2.0d0*ddr1(i)-ddr1(i+1)+d*(ddr1(i)+ddr1(i+1)
     3         +2.0d0*(dr1(i)-dr1(i+1)))))
        end if
        rhonudot=rhonu*adotoa*rhonudot/dlna
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function dtauda(a)
        implicit double precision (a-h,o-z)
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     &    grhonr
        common /qparm/ wdyn,ndyn
c
c  Compute expansion rate.
        if (a.lt.1.0d-9) then
           rhonu=1.0d0
           pnu=1.0d0/3.0d0
        else
           call nu1(a,rhonu,pnu)
        end if

c  8*pi*G*rho*a**4.
        if (ndyn.eq.0) then
                omegavdyn = omegav
        else
                omegavdyn = omegav * dynrho(a)
        endif

        grho2=grhom*(omegac+omegab)*a+grhog+grhor*annur+
     2  grhonr*annunr*rhonu+grhom*omegavdyn*a**4+grhom*omegak*a**2

#ifdef DIM

        grho2=(sqrt(grho2)+sqrt(omegav*grhom*a**4))**2
     1  +grhom*omegak*a**2
#endif

        dtauda=sqrt(3.0d0/grho2)

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function dsoundda(a)
        implicit double precision (a-h,o-z)
c
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     &                   grhonr
        double precision wdyn
        integer ndyn
        common /qparm/ wdyn,ndyn

c
c  Compute expansion rate.
        if (a.lt.1.0d-9) then
           rhonu=1.0d0
           pnu=1.0d0/3.0d0
        else
           call nu1(a,rhonu,pnu)
        end if
c  8*pi*G*rho*a**4.
        if (ndyn.eq.0) then
                omegavdyn = omegav
        else
                omegavdyn = omegav * dynrho(a)
        endif

        grho2=grhom*(omegac+omegab)*a+grhog+grhor*annur+
     2  grhonr*annunr*rhonu+grhom*omegavdyn*a**4+grhom*omegak*a**2

#ifdef DIM

        grho2=(sqrt(grho2)+sqrt(omegav*grhom*a**4))**2
     1  +grhom*omegak*a**2
#endif

        dtauda=sqrt(3.0d0/grho2)
c
        R=3.0d4*a*omegab*(h0/100.0d0)**2
        cs=1.0d0/sqrt(3.0d0*(1.0d0+R))
        dsoundda=dtauda*cs
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine reiopar(optdlss,zri,zristp,rif)

c This subroutine calculates the redshift of reionization
c from an optdlss and sets the reionization fraction rif=1.

      implicit double precision(a-h,o-z)

      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr

      zri=0.0d0
      rif=0.0d0
      if (optdlss.ne.0) then
         akthom=(2.3048d-9)*(1-yhe)*omegab*h0*h0
         rif=1.0d0

c Calculating the redshift of reionation the parameter
c needed by CMBFAST.
         na=1
         da=0.00001d0
         optd=0.0d0
 5       a=1-na*da
         optd=optd+da*rif*akthom*dtauda(a)/a**2
         if (optd.lt.optdlss) then
            na=na+1
            goto 5
         end if
         zri=1.0d0/a-1.0d0
         zristp=0.07d0*zri-1.0d0
         if (zristp.lt.0.0) zristp=0.0d0
c       TESTING
c        write(*,*)'zristp=',zristp
c       Problems in choosing zristp may lead to inaccuracies,
c       so it is better to always choose zristp=0, this leads
c       to a small increase in computation time. (Bug reported
c       by Andrew Liddle) July 2000.
         zristp=0.0

      end if
      return
      end

      subroutine recint(a,xe)
c calculates recfast xe
      implicit double precision (a-h,o-z)
      parameter (Nz0=10000)
      double precision zrec(Nz0),xrec(Nz0),dxrec(Nz0)
      common/recout/zrec,xrec,dxrec,Nz
      z=1.0/a-1.0
      if (z.ge.zrec(1)) then
       xe=xrec(1)
      else
       if (z.le.zrec(nz)) then
        xe=xrec(nz)
       else
        ilo=nz-z
        ihi=ilo+1
        az=z-int(z)
        bz=1.0d0-az
        xe=az*xrec(ilo)+bz*xrec(ihi)+
     1  ((az**3-az)*dxrec(ilo)+(bz**3-bz)*dxrec(ihi))/6.0d0
       endif
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine splder(y,dy,n)
c  Splder fits a cubic spline to y and returns the first derivatives at
c  the grid points in dy.  Dy is equivalent to a 4th-order Pade
c  difference formula for dy/di.
c
        implicit double precision (a-h,o-z)
        dimension f(40001),g(40001),y(n),dy(n)
        common /spl/ g
        save /spl/
c
        n1=n-1
        if (n1.gt.40000) write(*,*) 'Spline array overflow!!! n1=',
     2      n1,'>40000'
c  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
        f(1)=(-10.0d0*y(1)+15.0d0*y(2)-6.0d0*y(3)+y(4))/6.0d0
        f(n)=(10.0d0*y(n)-15.0d0*y(n1)+6.0d0*y(n-2)-y(n-3))/6.0d0
c  Solve the tridiagonal system
c  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
c  with dy(1)=f(1), dy(n)=f(n).
          do 10 i=2,n1
          f(i)=g(i)*(3.0d0*(y(i+1)-y(i-1))-f(i-1))
10      continue
        dy(n)=f(n)
          do 20 i=n1,1,-1
          dy(i)=f(i)-g(i)*dy(i+1)
20      continue
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine splini
c  Splini must be called before splder to initialize array g in common.
c
        implicit double precision (a-h,o-z)

        dimension g(40001)
        common /spl/ g
        save /spl/
c
        g(1)=0.0d0
          do 10 i=2,40001
          g(i)=1.0d0/(4.0d0-g(i-1))
10      continue
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
        parameter (MAXITER=40,MAXJ=5)
        implicit double precision (a-h,o-z)
        dimension g(MAXJ+1)
        double precision f
        external f
c
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
c  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      double precision yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100010)
      INTEGER i,k
      double precision p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        subroutine splint(y,z,n)
c  Splint integrates a cubic spline, providing the ouput value
c  z = integral from 1 to n of s(i)di, where s(i) is the spline fit
c  to y(i).
c
        implicit double precision (a-h,o-z)
        dimension y(n)
c
        n1=n-1
c  Cubic fit to dy/di at boundaries.
c       dy1=(-11.0d0*y(1)+18.0d0*y(2)-9.0d0*y(3)+2.0d0*y(4))/6.0d0
        dy1=0.0d0
        dyn=(11.0d0*y(n)-18.0d0*y(n1)
     1       +9.0d0*y(n-2)-2.0d0*y(n-3))/6.0d0
c
        z=0.5d0*(y(1)+y(n))+(dy1-dyn)/12.0d0
          do 10 i=2,n1
          z=z+y(i)
10      continue
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        SUBROUTINE output_power(ntf,amnu)
        include 'cmbfast.inc'
#ifdef MPIBUILD
        INCLUDE 'mpif.h'
        integer myid,sender,numprocs,ierr,status(MPI_STATUS_SIZE)
#endif
        integer ntf
        double precision amnu
        integer itf,iproc,j,icount
        integer ncount
        integer nmax
        parameter(nmax=1000)
        integer kmax
        parameter(kmax=1000)
        double precision datar(6*kmax)
        double precision k(kmax),tfc(kmax),tfb(kmax),tfg(kmax),
     >                  tfn(kmax),tfnm(kmax)
        integer indx(kmax)
        double precision power(nk0*6,ntfmax)
        integer istore(ntfmax)
        common/transfer_function/power,istore

#ifdef MPIBUILD
        CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )          !MPI
        CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
#else
        myid=0
        numprocs=1
#endif
        do itf = 1,ntf
#ifdef MPIBUILD
        CALL MPI_BARRIER_(MPI_COMM_WORLD,ierr)
#endif
        IF(MYID.ne.0) THEN
                if(istore(itf).gt.nmax) then
                        write(*,*) ' need to increase nmax'
#ifdef MPIBUILD
                        CALL MPI_ABORT_(MPI_COMM_WORLD,1)
#endif
                        stop
                endif
#ifdef MPIBUILD
                CALL MPI_SEND_(power(1,itf),6*nmax,
     $               MPI_DOUBLE_PRECISION,
     >               0,0,MPI_COMM_WORLD,IERR)
#endif
        ELSE
                do ncount = 1,istore(itf)
                        k(ncount) = power(6*(ncount-1)+1,itf)
                        tfc(ncount) = power(6*(ncount-1)+2,itf)
                        tfb(ncount) = power(6*(ncount-1)+3,itf)
                        tfg(ncount) = power(6*(ncount-1)+4,itf)
                        tfn(ncount) = power(6*(ncount-1)+5,itf)
                        tfnm(ncount) = power(6*(ncount-1)+6,itf)
                end do
            do iproc = 1,numprocs-1
#ifdef MPIBUILD
           CALL MPI_RECV_(datar,6*nmax,MPI_DOUBLE_PRECISION,
     >                          MPI_ANY_SOURCE,MPI_ANY_TAG,
     >                          MPI_COMM_WORLD,STATUS,IERR)
#endif
                        do j = 1,nmax
                                if(datar(6*(j-1)+1).ne.0.0) then
                                        ncount =ncount+1
                                else
                                        go to 20
                                endif
                                if(ncount.gt.kmax) then
                                   write(*,*) 'need to
     1                                  increase kmax'
#ifdef MPIBUILD
                                   CALL MPI_ABORT_(MPI_COMM_WORLD,1)
#endif
                                   STOP
                                endif
                                k(ncount) = datar(6*(j-1)+1)
                                tfc(ncount) = datar(6*(j-1)+2)
                                tfb(ncount) = datar(6*(j-1)+3)
                                tfg(ncount) = datar(6*(j-1)+4)
                                tfn(ncount) = datar(6*(j-1)+5)
                                tfnm(ncount) = datar(6*(j-1)+6)
                        end do !j
20                      continue
                        do j = 1,6*nmax
                                datar(j) = 0.0
                        end do
                end do !iproc
c
c  sort k's in increasing order
c
                call indexx(istore(itf),k,indx) !NR routine

                do icount=1,istore(itf)
                if(k(indx(icount)).ne.0.) then
                if (amnu.ne.0.0d0) then
                        write(12+itf,'(6E14.6)')k(indx(icount)),
     >                          tfc(indx(icount)),tfb(indx(icount)),
     >                          tfg(indx(icount)),tfn(indx(icount)),
     >                          tfnm(indx(icount))
                else
                        write(12+itf,'(5E14.6)')k(indx(icount)),
     >                          tfc(indx(icount)),tfb(indx(icount)),
     >                          tfg(indx(icount)),tfn(indx(icount))
                endif
                endif
                end do
        ENDIF
c       close(12+itf)
        END DO !itf

        RETURN
        END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine outtransf(n,y,curv,ik,itf)

        implicit double precision(a-h,o-z)

        dimension y(n)

        include 'cmbfast.inc'
        double precision power(nk0*6,ntfmax)
        integer istore(ntfmax)
        common/transfer_function/power,istore

        common /cosmoparm/ beta,beta2,amnu,lmaxg,lmaxnr,
     &                              lmaxnu,nqmax,iq0,iq1,iq2
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0,
     &                     tcmb,yhe,annur,annunr


        ak=sqrt(beta2-curv)
        ak2=ak*ak

        tfc=y(4)/ak2
        tfb=y(6)/ak2
        tfg=y(8)/ak2
        tfn=y(10+2*lmaxg)/ak2


        if (amnu.ne.0.0d0) then
           a=y(1)
           call nu1(a,rhonu,pnu)
           call nu2(a,drhonu,fnu,dpnu,shearnu,y(iq0),y(iq1),y(iq2))
           deltan=drhonu/rhonu
           tfnm=deltan/ak2
        endif

#ifdef MAX

c       FOR PROJECT WITH MAX
        omega0=omegac+omegab+omegan
        if (amnu.eq.0.0d0) then
           tf=(omegac*tfc+omegab*tfb)/omega0
!          write(12+itf,'(2E14.6)')beta,tf
        else
           tf=(omegac*tfc+omegan*tfnm+omegab*tfb)/omega0
!          write(12+itf,'(2E14.6)')beta,tf
        endif
#else

           power(6*istore(itf)+1,itf) =beta/(h0/100.d0)
           power(6*istore(itf)+2,itf) =tfc
           power(6*istore(itf)+3,itf) =tfb
           power(6*istore(itf)+4,itf) =tfg
           power(6*istore(itf)+5,itf) =tfn
        if (amnu.eq.0.0d0) then
c output k in units of Mpc/h

!          write(12+itf,'(5E14.6)')ak/(h0/100.d0),tfc,tfb,tfg,tfn

        else
!          write(12+itf,'(6E14.6)')ak/(h0/100.d0),tfc,tfb,tfg,
!     1       tfn,tfnm
        endif
        if (amnu.eq.0.0d0) then
                power(6*istore(itf)+6,itf) =0.0
        else
                power(6*istore(itf)+6,itf) =tfnm
        endif
        istore(itf) = istore(itf) +1

#endif

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE setuptransf(ntf)
        include 'cmbfast.inc'
        integer ntf
        integer i,itf
        double precision power(nk0*6,ntfmax)
        integer istore(ntfmax)
        common/transfer_function/power,istore
        if(ntf.gt.20) then
                write(*,*) ' need to increase dimensions in COMMON'
                write(*,*) ' block TRANSFER'
#ifdef MPIBUILD
                CALL MPI_FINALIZE_(ierr)
#endif
                stop
        endif
        do itf = 1,ntf
        istore(itf) = 0
        do i = 1,nk0*6
                power(i,itf)=0.0
        end do
        end do
        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine out(a,b,c,n,file)

c This subroutine just writes vectors a and b
c of lenght n to file

        double precision a(n),b(n),c(n)
        integer n,file,i

        do i=2,n
           write(file,*) i,a(i),b(i),c(i)
        end do

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCccc
C NUMERICAL RECIPES ROUTINE
C
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

c below are dark energy/quintessence routines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function dynrho(a)
        implicit double precision (a-h,o-z)
        include 'cmbfast.inc'

        common /qparm/ wdyn,ndyn
        double precision aw(nq0),wd(nq0),dwda(nq0)
        double precision wdpr(nq0),logrho(nq0)
        integer nq
        common /wvar/ aw,wd,wdpr,logrho,dwda,da,nq


       if(ndyn.ne.0) then

               if(ndyn.eq.1.or.ndyn.eq.3) then
                        dynrho = a**(-3.0d0*(1.0d0 + wdyn))

                else
                  klo=(a-aw(1))/da+1
                  khi=klo+1
                  dlna=log(a/aw(klo))
                  dynrho=logrho(klo)-3.0*((1+wd(klo)-
     1            dwda(klo)*aw(klo))*dlna+dwda(klo)*(a-aw(klo)))

                  dynrho = exp(dynrho)

                end if

       else
                  dynrho=1.0d0

       endif


        return
        end

        function wdyn_func(a)
        implicit double precision (a-h,o-z)
        include 'cmbfast.inc'

c       equation of state for dark energy
        common /qparm/ wdyn,ndyn
        double precision aw(nq0),wd(nq0),dwda(nq0)
        double precision wdpr(nq0),logrho(nq0)
        common /wvar/ aw,wd,wdpr,logrho,dwda,da,nq

c usual Omega_Lambda ndyn=0

       if(ndyn.ne.0) then

c  constant equation of state ndyn=1
c  tables with a,w,dw/da  ndyn=2

                if(ndyn.eq.1.or.ndyn.eq.3) then
                        wdyn_func = wdyn

                else if(ndyn.eq.2.or.ndyn.eq.4) then
                  klo=(a-aw(1))/da+1
                  khi=klo+1
                  as=(aw(khi)-a)/da
                  bs=1.0-as
                  wdyn_func=as*wd(klo)+bs*wd(khi)
c     2            +((as**3-as)*wdpr(klo)+
c     3            (bs**3-bs)*wdpr(khi))*(da**2)/6.d0
                 end if
       else

              wdyn_func=-1.0d0

       endif

        return
        end

        function dwda_func(a)
        implicit double precision (a-h,o-z)
        include 'cmbfast.inc'

c       equation of state for dark energy
        common /qparm/ wdyn,ndyn
        double precision aw(nq0),wd(nq0),dwda(nq0)
        double precision wdpr(nq0),logrho(nq0)
        common /wvar/ aw,wd,wdpr,logrho,dwda,da,nq

c usual Omega_Lambda ndyn=0

       if(ndyn.ne.0) then

c  constant equation of state ndyn=1
c  table with w(a)  ndyn=2

                if(ndyn.eq.1.or.ndyn.eq.3) then
                        dwda_func = 0.0

                else if(ndyn.eq.2.or.ndyn.eq.4) then
                  klo=(a-aw(1))/da+1
                  khi=klo+1
                  as=(aw(khi)-a)/da
                  bs=1.0-as
                  dwda_func=(wd(khi)-wd(klo))/da
c     2            +(-(3*as*as-1)*wdpr(klo)+
c     3            (3*bs*bs-1)*wdpr(khi))*da/6.d0
                 end if
       else
                 dwda_func = 0.0

       endif

        return
        end

        subroutine dyn_phi(a,hdot,grho,gpres,phi,psi,dphi,dpsi)

        implicit double precision (a-h,o-z)
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                            nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     &                 grhonr
        common /qparm/ wdyn,ndyn

        grhophi = dynrho(a) * omegav * grhom

        phid = sqrt(grhophi)*a

        source_term = -0.5*hdot*phid

        adot = sqrt(grho/3.0d0)

        addot = (grho - 3.0d0*gpres)/6.0d0

        coef1 = 2.0d0*adot

        if(ndyn.eq.1.or.ndyn.eq.3) then
           coef2 = -1.5d0*(1.0d0-wdyn)*(addot
     2         - adot*adot*(3.5d0+1.5d0*wdyn)) + ak2
           else
              wdyn=wdyn_func(a)
              dwda=dwda_func(a)
              if(wdyn.gt.-1.0) coef1=coef1+dwda*adot*a/(wdyn+1.0d0)

              coef2 = -1.5d0*(1.0d0-wdyn)*(addot -
     2          adot*adot*(3.5d0+1.5d0*wdyn)) + ak2 +
     3          3.0d0*(dwda*adot*a)*adot
        end if

           dphi = psi
           if (wdyn.gt.-1.0) then
            dpsi = -coef1*psi-coef2*phi+source_term
           else
            dpsi = 0.0d0
           endif

        return
        end

        subroutine dyn_nrg(a,grho,gpres,delphi,delpsi,
     1       dgrho,dgprs,dgtheta)

        implicit double precision (a-h,o-z)
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &                     ,tcmb,yhe,annur,annunr
        common /cosmoparm/ ak,ak2,amnu,lmax,lmaxnr,lmaxnu,
     &                            nqmax,iq0,iq1,iq2
        common /genparm/ grhom,grhog,grhor,adotrad,taurst,dtaurec,
     &          grhonr
        common /qparm/ wdyn,ndyn

        grhophi = dynrho(a) * omegav * grhom

        dphi = sqrt(grhophi)*a

        adot = sqrt(grho/3.0)

        if(ndyn .eq. 0) then
           dgrho = 0.0d0
           dgprs = 0.0d0
           dgtheta = 0.0d0
           return
        end if

        if(ndyn.eq.1.or.ndyn.eq.3) then

           dgrho = (1.0d0+wdyn)*dphi*(delpsi
     2          - 1.5d0*adot*(1.0d0-wdyn)*delphi)
           dgprs = (1.0d0+wdyn)*dphi*(delpsi
     2          + 1.5d0*adot*(1.0d0-wdyn)*delphi)
           dgtheta =  (1.0d0+wdyn)*ak2*dphi*delphi

           else
             w=wdyn_func(a)
             dwda=dwda_func(a)
             if (w.ge.-1.0) then
              dgrho = (1.0d0+w)*dphi*(delpsi -
     2          1.5d0*adot*(1.0d0-w)*delphi)
              dgprs = (1.0d0+w)*dphi*(delpsi +
     2          1.5d0*adot*(1.0d0-w)*delphi)
     3          + dphi*dwda*adot*a*delphi
              dgtheta =  (1.0d0+w)*ak2*dphi*delphi
             endif

           end if

        return
        end

        subroutine readtable
c reads w(a) and makes spline table; input must be linear in a
        implicit double precision(a-h,o-z)
        include 'cmbfast.inc'
        double precision aw(nq0),wd(nq0),dwda(nq0)
        double precision wdpr(nq0),logrho(nq0)
        integer nq
        common /wvar/ aw,wd,wdpr,logrho,dwda,da,nq

        d0hi=1.0d32
        d0lo=1.0d32
        rewind 99
        i=1
10      read(99,*,end=20)aw(i),wd(i)
        if (i.gt.1) then
c make density table by integrating up piecewise linear w(a)
         dlna=log(aw(i)/aw(i-1))
         dwda(i-1)=(wd(i)-wd(i-1))/(aw(i)-aw(i-1))
         logrho(i)=logrho(i-1)-3.0*((1+wd(i-1)-
     1   dwda(i-1)*aw(i-1))*dlna+dwda(i-1)*(aw(i)-aw(i-1)))
        else
         logrho(i)=0.0
        endif

        i=i+1
        goto 10
20      continue
        close(99)
        nq=i-1
        if (nq.gt.nq0) then 
         write(*,*)'You must increase nq0 in cmfast.inc to',nq
         stop
        endif
        do i=1,nq
         logrho(i)=logrho(i)-logrho(nq)
        enddo

c no spline to be consistent with logrho solution, which assumes w
c is piecewise linear in a
c        call spline(aw,wd,nq,d0lo,d0hi,wdpr)

        da=aw(2)-aw(1)

       return
       end
