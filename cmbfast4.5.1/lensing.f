C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (matiasz@sns.ias.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lensing(clts,cles,clbs,clcs,clkk,cltk)
      implicit double precision (a-h,o-z)

      include 'cmbfast.inc'

c      parameter (nnmax=10)

      double precision clts(l0max,nnmax)
      double precision cles(l0max,nnmax)
      double precision clbs(l0max,nnmax)
      double precision clcs(l0max,nnmax)
      double precision clkk(l0max,nnmax),cltk(l0max,nnmax)

      integer l(lmax),l0
      common /lvalues1/ l,l0,lmo
      save /lvalues1/

      common /initialps/ an(nnmax),alphans(nnmax),
     $     dalphansdlnk(nnmax),nn

      double precision s8(ntfmax,nnmax),anorm(ntfmax,nnmax)
      common /transfernorm/ s8,anorm

      common /ps/ xn,xnorm8,in
      common /idum/ ilin
      common /trlens/ntf
      common /lens/ lensflag

        double precision wdyn
        common /qparm/ wdyn,ndyn

c      double precision sigma,c2,integ
c      external sigma,c2,integ
      double precision sigma,c2
      external sigma,c2
      double precision bessj0,bessj1,bessj
      external bessj0,bessj1,bessj

c Decide if linear or non-linear
      ilin=0
      if (lensflag.eq.2) ilin=1

c Loop over power spectra
      do in=1,nn
         xn=an(in)
c Take the last transfer function as the one for today.
         xnorm8=anorm(ntf,in)
c Calculate epsilon and c2


         call epsilongen
c Lens the Cls


         call cllens(clts(1,in),cles(1,in),clbs(1,in)
     &        ,clcs(1,in),clkk(1,in),cltk(1,in),lmo)

      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine epsilongen
C
c     fluctuations in angular separation due to the potential fluctuations
c     from Limber's eqn. using nonlinear evolution. normalize to sigma8
c     units are Mpc
      implicit double precision (a-h,o-z)
      parameter (tol=1.0d-5,nx=1000,ntheta=50,hc=2.998d3)
      parameter (pi=3.1415926535898d0)
      parameter (kmax=5000,ntau=50,nknl=50,delta=1.0d0)
      dimension akk(kmax),tf(kmax),tfpr(kmax)
      dimension atrg(ntau,4)
      dimension xtau(nx),xa(nx),xpr(nx)
      dimension epsilon(ntheta),theta(ntheta)
      dimension c2oth(ntheta)
c      dimension epsprime(ntheta),c2othprime(ntheta)
      dimension delta2nl(nknl,ntau),alnknl(nknl,ntau),d2nlp(nknl,ntau)

        double precision wdyn
        common /qparm/ wdyn,ndyn

      common /tfun/ akk,tf,tfpr,k0
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr
      common /ps/ xn,xnorm8,in
      common /cosmo2/ curv,tau0,r0,omegam
      common /idum/ ilin
      common /trlens/ntf
      common /nlps/ delta2nl,alnknl,d2nlp
      common /par/atrg,tht,taur,rr
      common /outputeps1/ theta,epsilon,c2oth
      common /outputeps2/ nthetaout
      common /fitneut/apcb,fcb,yfsok2,fnu

      character*500 filetf
      common /filesepsgen/ filetf
      double precision sigma,c2
      double precision psk
      external sigma, c2
      external psl
      double precision dtauda
      external dtauda

      fourpi=4.0d0*pi
      nthetaout=ntheta

      omegam=omegac+omegab+omegan
      h=h0/100.0d0

      zeqp1=2.5d4*omegam*h*h*(2.7d0/tcmb)**4
c     EH (97) fitting formula for masive neutrino growth factor.
      fnu=omegan/omegam
      fcb=(omegac+omegab)/omegam
      if (fnu.gt.0.0d0) then
         apcb=0.25d0*(5.0d0-sqrt(1.0d0+24.0d0*fcb))
         aktoq=(2.7d0/tcmb)*(2.7d0/tcmb)/(omegam*h*h)*h
         yfsok2=17.2d0*fnu*(1.0d0+0.488d0
     &        *exp(-7.0d0*log(fnu)/6.0d0))
     &        *(dble(annunr)*aktoq/fnu)**2
      else
         apcb=1.0d0
         yfsok2=0.0d0
      end if


c     We will use transfer function output from CMBFAST and fitting
c     formulae for the growth factor. If more accuracy is needed
c     one can obtain the transfer function as a function of time
c     numerically from CMBFAST.

      rewind 12+ntf
      if (omegan.gt.0.0) then
         i=1
 30      read(12+ntf,*,end=20)akk(i),tcdm,tfb,tfg,tfn,tfnm
         akk(i)=log(akk(i))
         tf(i)=(omegac*tcdm+omegab*tfb+omegan*tfnm)/omegam
         tf(i)=log(tf(i))
c         if (akk(i).gt.1.0) goto 20
         i=i+1
         goto 30
 20      continue
      else
         i=1
 50      read(12+ntf,*,end=40)akk(i),tcdm,tfb,tfg,tfn
         akk(i)=log(akk(i))
         tf(i)=(omegac*tcdm+omegab*tfb)/omegam
         tf(i)=log(tf(i))
c         if (akk(i).gt.5.0) goto 40
         i=i+1

         goto 50
 40      continue
      end if

      k0=i-1
      akkhi=1.0d40
      akklo=1.0d40

      call spline(akk,tf,k0,akkhi,akklo,tfpr)


c     conformal time today
      a0=1.0d-8
      a1=1.0d0
      tau0=h*rombint(dtauda,a0,a1,tol)
      ar=1.0d0/1100.0d0
      taur=h*rombint(dtauda,a0,ar,tol)
      curv=-omegak/hc/hc
      alnkmin=akk(1)
      alnkmax=akk(k0)
      if (curv.ne.0.0d0) rcurv=1.0d0/sqrt(abs(curv))

c     set up (a,tau,r,g,n_eff) table for ntime points equally spaced
c     in time; use interpolation between a and tau

      xaold=1.0d-8
      xtauold=0.0d0

      do 15 i=1,nx
         xa(i)=1.0d0*i/nx
         xtau(i)=xtauold+h*rombint(dtauda,xaold,xa(i),tol)
         xaold=xa(i)
         xtauold=xtau(i)
 15   continue


      xaplo=1.0d0/dtauda(xa(1))
      xaphi=1.0d0/dtauda(1.0d0)
      call spline(xtau,xa,nx,xaplo,xaphi,xpr)

      do 45 i=1,ntau
         tau=dble(i)*(tau0-taur)/dble(ntau)+taur
         call splint_v2(xtau,xa,xpr,nx,tau,a)
         r=tau0-tau
         if (curv.ne.0.0d0) then
            if (curv.gt.0.) then
               r=sin(sqrt(curv)*r)/sqrt(curv)
            else
               r=sinh(sqrt(-curv)*r)/sqrt(-curv)
            endif
         endif
         if (r.lt.0.0d0) r=0.0d0
         x=1.0d0+omegam*(1.0d0/a-1)+omegav*(a*a-1.0d0)
         f=exp(4.0d0*log(omegam/a/x)/7.0d0)
c     Fitting formula for linear growth of perturbations.
c     Extra (1+zeq)=zeqp1 so that same normalization as EH 97 to use
c     their fitting formula for massive neutrino growth rate.
c     cancels in other cases.


         g=zeqp1*2.5d0*omegam/(x*f+1.5d0*omegam/a+1.0d0-omegam-omegav)

c Ma, Caldwell, Bode, Wang fitting formula

        if(wdyn.lt.0.0d0) then

        omnow = omegam/(omegam + omegav*(a**(-3.0d0*wdyn)))
        oqnow = 1.0d0 - omnow

        term1 = -(0.255d0*oqnow + 0.366d0*log(omnow))
        term2 = -(0.00266d0*oqnow - 0.07d0*log(omnow))/wdyn
        term3 = -(0.305d0*oqnow + 0.266d0*log(omnow))*wdyn

        g = g*( (-wdyn)**(term1+term2+term3) )

        endif

         atrg(i,1)=tau
         atrg(i,2)=a
         atrg(i,3)=r
         atrg(i,4)=g
 45   continue


      if (curv.eq.0.) then
         rr=tau0-taur
      else if (curv.gt.0.) then
         rr=sin(sqrt(curv)*(tau0-taur))/sqrt(curv)
      else
         rr=sinh(sqrt(-curv)*(tau0-taur))/sqrt(-curv)
      endif

c     If non linear calculation requested,
c     nonlinear mapping of power spectrum
      if (ilin.eq.1) then

c     only k>kcrit for which linear Delta^2>0.1  today.
c     find kcrit between kmin and kmax.

         alnklo=alnkmin
         alnkhi=alnkmax
         delta2crit=0.1d0
 60      if (alnkhi-alnklo.gt.0.01) then
            alnk=(alnkhi+alnklo)/2.0d0

            delta2=fourpi*xnorm8*exp(alnk*3.0d0)*psl(alnk)

            if (delta2.gt.delta2crit) then
               alnkhi=alnk
            else
               alnklo=alnk
            endif
            goto 60
         endif
         alnkcrit=alnk
         dlnk=(alnkmax-alnkcrit)/(nknl-1.0d0)
         g0=atrg(ntau,4)
         g0c=g0

c     TESTING
c      open(unit=89,file='deltatest.dat',status='unknown',
c     $     form='formatted')
c      open(unit=90,file='deltatest2.dat',status='unknown',
c     $     form='formatted')

         do 80 i=1,ntau

            g=atrg(i,4)
            gc=g
            a=atrg(i,2)

c            alnkl=alnkcrit
c            alnklow=alnkl-delta
c            t2low=log(psl(alnklow))
c            alnkhigh=alnkl+delta
c            t2high=log(psl(alnkhigh))
c            aneffold=(t2high-t2low)/2.0d0/delta

            do 70 j=1,nknl
               alnkl=alnkcrit+(j-1.0d0)*dlnk

c     EH (97) fitting formula for massive neutrinos
               if (fnu.gt.0.0d0) then
                  ak=exp(alnkl)
                  yfsp1=yfsok2*ak*ak+1.0d0
                  aux1=exp(0.7d0*log(g/yfsp1))
                  aux3=exp(0.7d0*log(fcb)/apcb)
                  aux1=aux1+aux3
                  aux2=apcb*log(aux1)/0.7d0
     &                 +(1.0d0-apcb)*log(g)
                  gc=exp(aux2)

                  aux1=exp(0.7d0*log(g0/yfsp1))
                  aux1=aux1+aux3
                  aux2=apcb*log(aux1)/0.7d0
     &                 +(1.0d0-apcb)*log(g0)
                  g0c=exp(aux2)
               end if

               delta2=fourpi*xnorm8*exp(3*alnkl)
     $         *psl(alnkl)*gc*gc/g0c/g0c

               alnklow=alnkl-delta-log(2.0d0)
               t2low=log(psl(alnklow))
               alnkhigh=alnkl+delta-log(2.0d0)
               t2high=log(psl(alnkhigh))
               aneff=(t2high-t2low)/2.0d0/delta
c               aneff=min(aneff,aneffold)
c               aneffold=aneff
c               if (i.eq.ntau) then
c                  write(90,'(8E15.5)')alnkl,log(psl(alnkl))
c     $                 ,alnklow,t2low,
c     $                 alnkhigh,t2high,aneff,aneff0
c               end if

               delta2nl(j,i)=fnl(delta2,g/a/zeqp1,aneff)
               alnknl(j,i)=log(1.0d0+delta2nl(j,i))/3.0d0+alnkl

c             write(35,'(4E12.4)')alnkl,delta2,alnknl(j,i),delta2nl(j,i)
c     use log(delta2) to facilitate interpolation
               delta2nl(j,i)=log(delta2nl(j,i))

 70         continue


            d2nlplo=1.d40
            d2nlphi=1.d40
            call spline(alnknl(1,i),delta2nl(1,i),nknl,d2nlplo,
     2           d2nlphi,d2nlp(1,i))
 80      continue
      end if

c     TESTING
c      close(89)
c      close(90)
c
c      open(unit=89,file='deltanl.dat',status='unknown',
c     $     form='formatted')
c      do j=1,nknl
c         write(89,'(11E15.5)')alnknl(j,ntau),delta2nl(j,ntau),
c     $      alnknl(j,ntau-1),delta2nl(j,ntau-1),delta2nl(j,ntau-2),
c     $      delta2nl(j,ntau-5),delta2nl(j,ntau-7),delta2nl(j,ntau-9),
c     $      delta2nl(j,ntau-11),delta2nl(j,ntau-13),delta2nl(j,ntau-15)
c      end do
c      close(89)


c     k integration
      theta0=pi/180.d0/600.d0
      theta1=0.8d0
      dlntheta=(log(theta1)-log(theta0))/(ntheta-1)
c     TESTING
c      write(*,*)'ntheta',ntheta,alnkmin,alnkmax
c      open(unit=95,file='sigma.dat',status='unknown',
c     $     form='formatted')
c      open(unit=661,file='sigmadit.dat',status='unknown',
c     $     form='formatted')
c      ntest=1000
c      tht=theta0
c      do j=1,ntest
c         write(*,*)'j sigma=',j
c         alnktest=alnkmin+dble(j-1)*(alnkmax-alnkmin)/dble(ntest-1)
c         sigmatest=sigma(alnktest)
c         write(95,'(2E15.5,1I6)')alnktest,sigmatest,j
c      end do
c      close(661)
c      close(95)


      do 130 itheta=1,ntheta

         theta(itheta)=exp(log(theta0)+dlntheta*(itheta-1))
         tht=theta(itheta)
         sigth=sqrt(fourpi*rombint(sigma,alnkmin,alnkmax,tol))
         xkappa=1/tht
         sigma2kappa=xkappa*xkappa*fourpi*psk(xkappa)
c         write(33,*)xkappa,sigma2kappa/4.0
         epsilon(itheta)=sigth/theta(itheta)
c         c2th=sqrt(fourpi*rombint(c2,alnkmin,alnkmax,tol))
         c2th=(fourpi*rombint(c2,alnkmin,alnkmax,tol))
         c2oth(itheta)=c2th/theta(itheta)/theta(itheta)
 130  continue


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       function fnl(x,g,aneff)
       implicit double precision (a-h,o-z)
c nonlinear mapping from Peacock and Dodds
c quick fix for WDM, really not appropriate
       if (aneff.gt.-2.99d0) then
         dum=log(1.0d0+aneff/3.0d0)
       else
         dum=-5.7
       endif
       a=0.482d0*exp(-0.947d0*dum)
       b=0.226d0*exp(-1.778d0*dum)
       alpha=3.310d0*exp(-0.244d0*dum)
       beta=0.862d0*exp(-0.287d0*dum)
       v=11.55d0*exp(-0.423d0*dum)
       dum1=1.0d0+b*beta*x+exp(alpha*beta*log(a*x))
       dum2=1.0d0+exp(beta*(alpha*log(a*x)+3.0d0*log(g)-
     2      log(v)-0.5d0*log(x)))
       fnl=x*exp(log(dum1/dum2)/beta)
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function psk(xkappa)
      implicit double precision (a-h,o-z)
      parameter (ntau=50,hc=2.998d3,nknl=50)
      dimension dum1(ntau)
      dimension atrg(ntau,4)
      dimension delta2nl(nknl,ntau),alnknl(nknl,ntau),d2nlp(nknl,ntau)
      common /par/atrg,theta,taur,rr
      common /ps/ xn,xnorm8,in
      common /cosmo2/ curv,tau0,r0,omegam
      common /idum/ ilin
      common /nlps/ delta2nl,alnknl,d2nlp
      common /fitneut/apcb,fcb,yfsok2,fnu
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr

c     computes integrand of sigma as a function of theta
c     do the time integration at a fixed number of points
c     using preevaluated nonlinear power spectrum; interpolate
c     to get the ps at a given k

      h=h0/100.0d0

      fourpi=4.0d0*3.1415926535898d0
      g0=atrg(ntau,4)
      curv=-omegak/hc/hc
      if (curv.ne.0.0d0) rcurv=1.0d0/sqrt(abs(curv))


c EH (97) fitting formula for massive neutrinos

      if (fnu.gt.0.0d0) then
         yfsp1=yfsok2*ak*ak+1.0d0
         aux1=exp(0.7d0*log(g0/yfsp1))
         aux1=aux1+exp(0.7d0*log(fcb)/apcb)
         aux2=apcb*log(aux1)/0.7d0
     &        +(1.0d0-apcb)*log(g0)
         g0=exp(aux2)
      end if

      do 10 i=1,ntau-1

         a=atrg(i,2)
         r=atrg(i,3)
         g=atrg(i,4)
         ak=xkappa/r
         alnk=log(ak)


c EH (97) fitting formula for massive neutrinos
         if (fnu.gt.0.0d0) then
            yfsp1=yfsok2*ak*ak+1.0d0
            aux1=exp(0.7d0*log(g/yfsp1))
            aux1=aux1+exp(0.7d0*log(fcb)/apcb)
            aux2=apcb*log(aux1)/0.7d0
     &           +(1.0d0-apcb)*log(g)
            g=exp(aux2)
         end if
         x=ak*r*theta

         if (ilin.ne.1) then

            ps=9.0d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2           *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam
         else

c     linear theory next 2 lines
            if (alnk.lt.alnknl(1,i).or.alnk.gt.alnknl(nknl,i)) then
               ps=9.d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2              *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam

            else
c     nonlinear theory next 4 lines

               call splint_v2(alnknl(1,i),delta2nl(1,i),d2nlp(1,i),
     2              nknl,alnk,delta2)
               ps=9.0d0/4.0d0/fourpi*exp(delta2)*exp(-4.0d0*alnk)
     2              /hc/hc/hc/hc/a/a*omegam*omegam
            endif
         endif

         w=sqrt(1.0d0-curv*r*r)-sqrt(1.0d0-curv*rr*rr)*r/rr

         dum1(i)=fourpi*ps*ak*w*w
 10   continue

      dum1(ntau) = 0

      dtau=(tau0-taur)/ntau

c     integrate in equal time intervals

      call integ(dum1,psk,ntau)
      psk=psk*dtau
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function sigma(alnk)
      implicit double precision (a-h,o-z)
      parameter (ntau=50,hc=2.998d3,nknl=50)
      dimension dum1(ntau)
      dimension atrg(ntau,4)
      dimension delta2nl(nknl,ntau),alnknl(nknl,ntau),d2nlp(nknl,ntau)
      common /par/atrg,theta,taur,rr
      common /ps/ xn,xnorm8,in
      common /cosmo2/ curv,tau0,r0,omegam
      common /idum/ ilin
      common /nlps/ delta2nl,alnknl,d2nlp
      common /fitneut/apcb,fcb,yfsok2,fnu
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr

c     computes integrand of sigma as a function of theta
c     do the time integration at a fixed number of points
c     using preevaluated nonlinear power spectrum; interpolate
c     to get the ps at a given k

      h=h0/100.0d0

      fourpi=4.0d0*3.1415926535898d0
      ak=exp(alnk)
      g0=atrg(ntau,4)
      curv=-omegak/hc/hc
      if (curv.ne.0.0d0) rcurv=1.0d0/sqrt(abs(curv))


c EH (97) fitting formula for massive neutrinos

      if (fnu.gt.0.0d0) then
         yfsp1=yfsok2*ak*ak+1.0d0
         aux1=exp(0.7d0*log(g0/yfsp1))
         aux1=aux1+exp(0.7d0*log(fcb)/apcb)
         aux2=apcb*log(aux1)/0.7d0
     &        +(1.0d0-apcb)*log(g0)
         g0=exp(aux2)
      end if

      do 10 i=1,ntau

         a=atrg(i,2)
         r=atrg(i,3)
         g=atrg(i,4)

c EH (97) fitting formula for massive neutrinos
         if (fnu.gt.0.0d0) then
            yfsp1=yfsok2*ak*ak+1.0d0
            aux1=exp(0.7d0*log(g/yfsp1))
            aux1=aux1+exp(0.7d0*log(fcb)/apcb)
            aux2=apcb*log(aux1)/0.7d0
     &           +(1.0d0-apcb)*log(g)
            g=exp(aux2)
         end if
         x=ak*r*theta

         if (ilin.ne.1) then

            ps=9.0d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2           *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam
         else

c     linear theory next 2 lines
            if (alnk.lt.alnknl(1,i)) then
               ps=9.d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2              *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam

            else
c     nonlinear theory next 4 lines


               call splint_v2(alnknl(1,i),delta2nl(1,i),d2nlp(1,i),
     2              nknl,alnk,delta2)
               ps=9.0d0/4.0d0/fourpi*exp(delta2)*exp(-4.0d0*alnk)
     2              /hc/hc/hc/hc/a/a*omegam*omegam
            endif
         endif

         if (x.lt.0.1) then
            aj0=1.0d0
            daj0=x*x/4.0d0
         else
            aj0=bessj0(x)
            daj0=1.0d0-aj0
         endif
         w=sqrt(1.0d0-curv*r*r)-sqrt(1.0d0-curv*rr*rr)*r/rr
c     TESTING
c         write(661,'(5E15.5,1I6)')ak,ps,w,daj0,alnk,i

         dum1(i)=fourpi*ps*ak*w*w*daj0
 10   continue

      dtau=(tau0-taur)/ntau

c     integrate in equal time intervals

      call integ(dum1,sigma,ntau)
      sigma=sigma*dtau
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function c2(alnk)
      implicit double precision (a-h,o-z)
      parameter (ntau=50,hc=2.998d3,nknl=50)
      dimension dum1(ntau)
      dimension atrg(ntau,4)
      dimension delta2nl(nknl,ntau),alnknl(nknl,ntau),d2nlp(nknl,ntau)
      common /par/atrg,theta,taur,rr
      common /ps/ xn,xnorm8,in
      common /cosmo2/ curv,tau0,r0,omegam
      common /idum/ ilin
      common /nlps/ delta2nl,alnknl,d2nlp
      common /fitneut/apcb,fcb,yfsok2,fnu
      common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr

c     computes integrand of sigma as a function of theta
c     do the time integration at a fixed number of points
c     using preevaluated nonlinear power spectrum; interpolate
c     to get the ps at a given k

      h=h0/100.0d0

      fourpi=4.0d0*3.1415926535898d0
      ak=exp(alnk)
      g0=atrg(ntau,4)

      curv=-omegak/hc/hc
      if (curv.ne.0.0d0) rcurv=1.0d0/sqrt(abs(curv))

c EH (97) fitting formula for massive neutrinos

      if (fnu.gt.0.0d0) then
         yfsp1=yfsok2*ak*ak+1.0d0
         aux1=exp(0.7d0*log(g0/yfsp1))
         aux1=aux1+exp(0.7d0*log(fcb)/apcb)
         aux2=apcb*log(aux1)/0.7d0
     &        +(1.0d0-apcb)*log(g0)
         g0=exp(aux2)
      end if

      do 10 i=1,ntau

         a=atrg(i,2)
         r=atrg(i,3)
         g=atrg(i,4)

c EH (97) fitting formula for massive neutrinos
         if (fnu.gt.0.0d0) then
            yfsp1=yfsok2*ak*ak+1.0d0
            aux1=exp(0.7d0*log(g/yfsp1))
            aux1=aux1+exp(0.7d0*log(fcb)/apcb)
            aux2=apcb*log(aux1)/0.7d0
     &           +(1.0d0-apcb)*log(g)
            g=exp(aux2)
         end if

         x=ak*r*theta

         if (ilin.ne.1) then

            ps=9.0d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2           *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam
         else

c     linear theory next 2 lines
            if (alnk.lt.alnknl(1,i)) then
               ps=9.d0/4.0d0*xnorm8*psl(alnk)/exp(alnk)
     2              *g*g/g0/g0/hc/hc/hc/hc/a/a*omegam*omegam

            else
c     nonlinear theory next 4 lines


               call splint_v2(alnknl(1,i),delta2nl(1,i),d2nlp(1,i),
     2              nknl,alnk,delta2)
               ps=9.0d0/4.0d0/fourpi*exp(delta2)*exp(-4.0d0*alnk)
     2              /hc/hc/hc/hc/a/a*omegam*omegam
            endif
         endif

         if (x.lt.0.1) then
            aj2=x*x/8.0d0
         else
            aj2=bessj(2,x)
         endif
         w=sqrt(1.0d0-curv*r*r)-sqrt(1.0d0-curv*rr*rr)*r/rr
         dum1(i)=fourpi*ps*ak*w*w*aj2
 10   continue

      dtau=(tau0-taur)/ntau

c     integrate in equal time intervals

      call integ(dum1,c2,ntau)
      c2=c2*dtau

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function psl(alnk)
c  t evaluates the transfer function at wavenumber ak.
c
        implicit double precision (a-h,o-z)
        parameter (kmax=5000)
        dimension akk(kmax),tf(kmax),tfpr(kmax)
        common /tfun/ akk,tf,tfpr,k0
        common /ps/ xn,xnorm8,in
        common /lingerinc/ omegab,omegac,omegav,omegan,omegak,h0
     &     ,tcmb,yhe,annur,annunr
        h=h0/100.0d0
         if (alnk.lt.akk(1)) then
          t2=1.0d0
         else
          call splint_v2(akk,tf,tfpr,k0,alnk,t2)
          t2=exp(2.0d0*t2)
         endif
         ak=exp(alnk)*h
         if (abs(omegak).le.1.0d-3) then
             call powersflat(ak,in,apowers)
          else
             call powersopen(ak,in,apowers)
          end if
          psl=apowers*ak/h*t2
         return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cllens(clt,cle,clb,clc,clkk,cltk,lmaxin)
c     calculating the multipoles C_l after lensing using sigma(theta)
      implicit double precision (a-h,o-z)

      include 'cmbfast.inc'

      parameter (pi=3.141592654d0)
      double precision clt(l0max),cle(l0max),clb(l0max),clc(l0max)
      double precision cltl(l0max),clel(l0max),clbl(l0max),clcl(l0max)
      double precision clkk(l0max),cltk(l0max)
c      double precision alout(l0max)

      parameter (nth=50)
      dimension epsilon(nth),theta(nth),epsprime(nth)
      dimension c2oth(nth),c2othprime(nth)
      common /outputeps1/ theta,epsilon,c2oth
      common /outputeps2/ ntheta

c      parameter (nx=26000)
      parameter (nx=7*l0max+2)

      double precision aj0(nx),daj0(nx)
      double precision aj2(nx),daj2(nx)
      double precision aj4(nx),daj4(nx)
      double precision aj6(nx),daj6(nx)

      character*500 filejlens
      common /filej/ filejlens


c     Initializing Cl arrays
      do l=1,l0max
         cltl(l)=0.0d0
         clel(l)=0.0d0
         clbl(l)=0.0d0
         clcl(l)=0.0d0
      end do

c     Getting Cls
      do l=2,lmaxin
         clt(l)=clt(l)/dble(l)/(dble(l)+1.0d0)
         cle(l)=cle(l)/dble(l)/(dble(l)+1.0d0)
         clc(l)=clc(l)/dble(l)/(dble(l)+1.0d0)
      end do


c     Get interpolation tables for epsilon and c2.
      ith=ntheta
      theta(ith)=pi
      epsilon(ith)=0.0d0
      c2oth(ith)=0.0d0
      epslo=1.d40
      epshi=1.d40
      call spline(theta,epsilon,ith,epslo,
     2     epshi,epsprime)
      call spline(theta,c2oth,ith,epslo,
     2     epshi,c2othprime)


c     Read array with bessel function

      open(unit=20,file=filejlens
     <     ,status='old',form='unformatted')
      rewind 20
      read(20)nxfile
      read(20)lmaximfile
      read(20)nt
      read(20)dx
      read(20)aj0
      read(20)daj0
      read(20)aj2
      read(20)daj2
      read(20)aj4
      read(20)daj4
      read(20)aj6
      read(20)daj6
      close(20)

      if (nxfile.ne.nx) then
         write(*,*)'jlens.dat and lensing.f have inconsistent'
         write(*,*)'l0max, you will have to rerun jlens.dat'
         stop
      end if
      xmaxlens=pi*dble(lmaxin)/10.0
      if (lmaxin.gt.lmaximfile) then
         write(*,*)'jlens.dat was was run with a lower lmax'
         write(*,*)' than needed. You will have to rerun it'
      end if

c     No need to integrate to pi because we are calculating
c     the deviations due to lensing which are very small on large
c     angular scales

      l0=2*lmaxin
      m=l0/10
      dtheta=pi/dble(l0)
      do 40 i=1,m

         theta1=dble(i)*dtheta

         call splint_v2(theta,epsilon,epsprime,ith,theta1,eps1)
         sig1=(eps1*theta1)**2

         call splint_v2(theta,c2oth,c2othprime,ith,theta1,c2gl)
         c2gl=c2gl*theta1*theta1

c     get the correlation functions at theta
         cth1=0.0d0
         cth2=0.0d0
         cth3=0.0d0
         cth4=0.0d0

         do 60 l1=2,lmaxin
            x1=dble(l1)*theta1

c     Interpolate Bessel function

            i1lo=int(x1/dx)+1
            i1hi=i1lo+1
            a=dble(i1lo)-1.0d0-x1/dx
            b=1.0d0-a
            a31=(a**3-a)*dx*dx/6.0d0
            b31=(b**3-b)*dx*dx/6.0d0

            b0=a*aj0(i1lo)+b*aj0(i1hi)+a31*daj0(i1lo)
     &           +b31*daj0(i1hi)

            b2=a*aj2(i1lo)+b*aj2(i1hi)+a31*daj2(i1lo)
     &           +b31*daj2(i1hi)

            b4=a*aj4(i1lo)+b*aj4(i1hi)+a31*daj4(i1lo)
     &           +b31*daj4(i1hi)

            b6=a*aj6(i1lo)+b*aj6(i1hi)+a31*daj6(i1lo)
     &           +b31*daj6(i1hi)

            aux1=exp(-0.5d0*sig1*dble(l1*l1))
            aux2=0.5d0*dble(l1*l1)*c2gl*aux1
            aux1=x1*(aux1-1.0d0)
            aux2=x1*aux2

c     cth4 is - the correlation function, when we go back to
c     fourier space that sign is compensated.

            cth1=cth1+(aux1*b0+aux2*b2)*clt(l1)
            cth2=cth2+(aux1*b0+aux2*b2)*cle(l1)
            cth3=cth3+(aux1*b4+aux2*0.5d0*(b2+b6))*cle(l1)
            cth4=cth4+(aux1*b2+aux2*0.5d0*(b0+b4))*clc(l1)


 60      continue

         do 70 l=2,lmaxin

            x2=dble(l)*theta1

            i1lo=int(x2/dx)+1
            i1hi=i1lo+1
            a=dble(i1lo)-1.0d0-x2/dx
            b=1.0d0-a
            a31=(a**3-a)*dx*dx/6.0d0
            b31=(b**3-b)*dx*dx/6.0d0

            b0=a*aj0(i1lo)+b*aj0(i1hi)+a31*daj0(i1lo)
     &           +b31*daj0(i1hi)

            b2=a*aj2(i1lo)+b*aj2(i1hi)+a31*daj2(i1lo)
     &           +b31*daj2(i1hi)

            b4=a*aj4(i1lo)+b*aj4(i1hi)+a31*daj4(i1lo)
     &           +b31*daj4(i1hi)

            cltl(l)=cltl(l)+b0*cth1
            clel(l)=clel(l)+0.5d0*(b0*cth2+b4*cth3)
            clbl(l)=clbl(l)+0.5d0*(b0*cth2-b4*cth3)
            clcl(l)=clcl(l)+b2*cth4

 70      continue
 40   continue

c     Last 300 ls are not accurate because of convolution.

      do 80 l=2,lmaxin-300

         aux1=dble((l+1)*l)
         aux2=aux1*dtheta

         clt(l)=cltl(l)*aux2+clt(l)*aux1
         cle(l)=clel(l)*aux2+cle(l)*aux1
         clb(l)=clbl(l)*aux2
         clc(l)=clcl(l)*aux2+clc(l)*aux1

         write(10,'(1I5,6E13.4)')l,clt(l)
     &        ,cle(l),clb(l),clc(l),clkk(l),cltk(l)


 80   continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine integ(y,z,n)
c  integ integrates using a n-order polynomial
c
        implicit double precision (a-h,o-z)
        dimension y(n)
c include tau=0 point with value 0
        z=y(n)*3.0d0/8.0d0+(y(1)+y(n-1))*7.0d0/6.0d0
     >                  +(y(2)+y(n-2))*23.0d0/24.0d0
        do 10 i=3,n-3
         z=z+y(i)
10      continue
        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION bessj0(x)
      DOUBLE PRECISION bessj0,x
      DOUBLE PRECISION ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.d0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION BESSJ(N,X)
      double precision bessj,x,bessj0,bessj1
      PARAMETER (IACC=40,BIGNO=1.E10,BIGNI=1.E-10)
      IF(N.LT.2)STOP 'bad argument N in BESSJ'
      TOX=2./real(X)
      IF(X.GT.FLOAT(N))THEN
        BJM=real(BESSJ0(X))
        BJ=real(BESSJ1(X))
        DO 11 J=1,N-1
          BJP=J*TOX*BJ-BJM
          BJM=BJ
          BJ=BJP
11      CONTINUE
        BESSJ=dble(BJ)
      ELSE
        M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
        BESSJ=0.0d0
        JSUM=0
        SUM=0.
        BJP=0.
        BJ=1.
        DO 12 J=M,1,-1
          BJM=J*TOX*BJ-BJP
          BJP=BJ
          BJ=BJM
          IF(ABS(BJ).GT.BIGNO)THEN
            BJ=BJ*BIGNI
            BJP=BJP*BIGNI
            BESSJ=BESSJ*dble(BIGNI)
            SUM=SUM*BIGNI
          ENDIF
          IF(JSUM.NE.0)SUM=SUM+BJ
          JSUM=1-JSUM
          IF(J.EQ.N)BESSJ=dble(BJP)
12      CONTINUE
        SUM=2.*SUM-BJ
        BESSJ=BESSJ/dble(SUM)
      ENDIF
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION BESSJ1(X)
      double precision bessj1,x,ax,z,xx
      DOUBLE PRECISION Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,
     *    R6,S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        rex=real(x)
        Z=8.0d0/AX
        Y=Z**2
        XX=AX-2.356194491d0
        BESSJ1=SQRT(.636619772d0/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *dble(SIGN(1.,rex))
      ENDIF
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE splint_v2(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.d0
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
