      program ujlgen
C     This code finds the initial conditions
c     for the ujl integration. It finds for each
c     beta and l the value of chi where ujl*l=10^{-6}.
c     This choice is motivated by the WKB approximation
c     for these functions as their amplitude scales with
c     l. This choice ensures that the initial condition
c     is in the regime where the potential
c     is much larger than the energy, in the wkb
c     approximation. In this new version, that includes
c     closed universes.
      implicit double precision(a-h,o-z)

      include 'cmbfast.inc'

      parameter (nbmax=1000)
      double precision x1,x2,yo1,yo2
      double precision xo,dx
      double precision y1,y2
      double precision rtsafe

      double precision ay2o(nbmax),axo(nbmax)
      double precision ay2c(nbmax),axc(nbmax)
      double precision abeta(nbmax)

      integer l(lmax)
      logical intflg
      character*500 filename,arg
      common /param/ beta,ll,K
      common /out/ yout1,yout2

      common /lvalues1/ l,l0,lmo
      save /lvalues1/

      external rtsafe
      external fujl

c     TIMING
c     Timing variables for testing purposes
c     Should be comented in released versions
c     real etime,actual,timeprev,tarray(2)
c     external etime

#ifdef IARGC
      if (iargc() .eq. 0) then
         intflg = .true.
      else
         intflg = .false.
         if (iargc() .ne. 3) then
            print *, 'Usage: ujlgen [lmax] [kmax] [output filename]'
            stop 1
         end if
         call getarg(1,arg)
         read (arg,*) lmoin
         call getarg(2,arg)
         read (arg,*) kmax0
         call getarg(3,filename)
      endif
#else
      intflg = .true.
#endif

      if (intflg) then
         write(*,*)'Maximum value of l,keta: (1500,3000)'
         read(*,*)lmoin,kmax0
      endif


      if ((lmoin+300).gt.l0max) then
         write(*,*)'Sorry there is a maximum l of ', l0max
         write(*,*)'The code needs to calculate up to',lmoin+300
         write(*,*)'300 more than what you asked.'
         write(*,*)'Please increase l0max in '
         write(*,*)'jlgen.f, ujlgen.f,'
         write(*,*)'jlens.f,ccmbflat.f and cmbopen.f'
         write(*,*)'in all places where it appears'
         stop
      end if

      call initlval(lmoin)

c     This file will contain both initial conditions
c     for open and closed models.
      if (intflg) then
         write(*,*) 'Enter output filename'
         read(*,'(a500)') filename
      endif
      open(unit=9,file=filename
     <     ,status='unknown',form='unformatted')
      rewind 9

      betamax=1.0d4*kmax0
      dlnbeta=20.0d0
      betamin=1.0d-4
      Npoint=500
      write(9)lmo
      write(9)betamax
      write(9)betamin
      write(9)dlnbeta
      write(9)Npoint

c     For safety purposes check. We will use this to
c     check if the lvalues are the same in file and code.
        write(9)l0
        do j=1,l0
           write(9)l(j)
        end do


C     Setting BETAS FOR OPEN MODELS
      ntot=int(log(betamax/betamin)*dlnbeta)
      do i=1,ntot
         abeta(i)=betamin*exp(dble(i-1)/dlnbeta)
      end do
      do il=1,l0

c     TIMING
c     timeprev=etime(tarray)

         ll=l(il)
         tol=max(1.0d-4/dble(ll),1.0d-8)

cxxxc     TESTING
cxxx         write(*,*)il,ll

C     FIRST OPEN MODELS
         K=-1

         do i=1,ntot

            beta=abeta(i)

            if (i.eq.1) then

C     FIRST beta for this l, use simple formula
C     to guess starting value

               alzp=19.736d0
               if(ll.eq.2)alzp=1258.0d0
               if(ll.eq.3)alzp=199.52d0
               if(ll.eq.4)alzp=79.43d0
               if(ll.eq.5)alzp=50.11d0
               if(ll.eq.6)alzp=35.48d0
               if(ll.eq.7)alzp=29.52d0

               wstart=beta/alzp
               xstart=wstart*sqrt(dble(ll*(ll+1)))
               x1=log(xstart)
               dx=0.5d0
c

            else

c     Use solution of previous beta as starting value
c     also use previous interval to guess some dx.
               x1=xo
               dx=abs((x1-x2)/(yo1-yo2))
            end if

c     TESTING. Output starting value and dx for one l
c     as a function of beta.
c            if (ll.eq.2) then
c               write(83,'(3E15.5)')beta,x1,dx
c            end if

C     Bracket Root
            yo2=0
            x2=x1
            call braroot(x1,x2,yo1,yo2,dx,fujl)

C     Now Root is between x1 and x2
C     Use rtsafe to find it.

            xo=rtsafe(fujl,x2,x1,tol)
            call fujl(xo,y1,y2,9)
            axo(i)=dble(xo)
            ay2o(i)=yout2

c     TESTING. Output root, ujl and derivative for one l
c     as a function of beta.
c            if (ll.eq.2) then
c               write(82,'(4E15.5)')beta,xo,yout1,yout2
c            end if

         end do


C     NOW CLOSED MODELS
         K=1

c     I will find the roots in reverse order in beta
c     for closed models. This is convenient because
c     at high beta the root for closed and open
c     expressed in terms of w coincide so I know
c     where to start.

c     The list of betas for closed models will be
c     l+1,l+2,...,l+Npoint,beta(ist+1),beta(ist+2),....,beta(ntot)
c     where ist is given by

         aist=1.0d0+dlnbeta*log(dble(ll+Npoint+1)/betamin)
         ist=int(aist)

c     So that the total numbers of betas for this l are
         ntotc=ntot-ist+Npoint

         do i=1,ntotc
            if (i.le.(ntotc-Npoint)) then
               beta=abeta(ntot-i+1)
               beta=dble(int(beta))
            else
               beta=dble(ll+1+ntotc-i)
            end if

c     Use solution of previous beta as starting value
c     also use previous interval to guess some dx.
c     For the first beta this is the solution coming
c     from open models.
            dx=abs((x1-x2)/(yo1-yo2))
            x1=xo

c     TESTING. Output starting value and dx for one l
c     as a function of beta.
c            if (ll.eq.2) then
c               write(81,'(3E15.5)')beta,x1,dx
c            end if

C     Bracket Root
            yo2=0
            x2=x1
            call braroot(x1,x2,yo1,yo2,dx,fujl)


C     Now Root is between x1 and x2
C     Use rtsafe to find it.

            xo=rtsafe(fujl,x2,x1,tol)

            call fujl(xo,y1,y2,9)

            axc(i)=dble(xo)
            ay2c(i)=yout2


c     TESTING. Output root, ujl and derivative for one l
c     as a function of beta.
c            if (ll.eq.10) then
c               write(80,'(4E15.5)')beta
c     $              ,exp(xo)/sqrt(dble(ll*(ll+1)))
c     $              ,yout1,yout2
c            end if
c            if (ll.eq.1600) then
c               write(81,'(4E15.5)')beta
c     $              ,exp(xo)/sqrt(dble(ll*(ll+1)))
c     $              ,yout1,yout2
c            end if
c
         end do

c     Output

c     Open models
         do i=1,ntot
            write(9)axo(i)
            write(9)ay2o(i)
         end do
c     Closed  models
         do i=1,ntotc
            write(9)axc(i)
            write(9)ay2c(i)
         end do

c     TIMING
c     actual=etime(tarray)
c     write(51,*)ll,actual-timeprev

      end do
      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine initlval(lmoin)

c This subroutines initializes l arrays.

        implicit double precision(a-h,o-z)

        include 'cmbfast.inc'

        integer l(lmax),l0

        common /lvalues1/ l,l0,lmo
        save /lvalues1/

c 300 more ls will always be calculated, this is necessary
c for lensing.

        lmo=lmoin+300

        include 'lvalues.inc'

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine braroot(x1,x2,yo1,yo2,dx,fujl)

c     After this routine is called the root should
c     be inside (x1,x2). yo1 and yo2 are the corresponding
c     values of fujl at x1 and x2. dx is the size of the
c     step used to bracket solution. When called x1 is some
c     starting value near the root, yo2 is set to 0.
c     Other things should be irrelvant.

      implicit double precision(a-h,o-z)
      external fujl
      common /param/ beta,ll,K

 50         continue
            call fujl(x1,yo1,yo3,5)
            if ((yo1*yo2).lt.0) goto 60

            if (yo1.lt.0) then
               x2=x1
               yo2=yo1
               x1=x1+dx
               goto 50
            end if

            if (yo1.gt.0) then
               x2=x1
               yo2=yo1
               x1=x1-dx
               goto 50
            end if
 60         continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION rtsafe(funcd,x1,x2,xacc)
      INTEGER MAXIT
      double precision rtsafe,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (MAXIT=1000)
      INTEGER j
      double precision df,dx,dxold,f,fh,fl,temp,xh,xl
      call funcd(x1,fl,df,10)
      call funcd(x2,fh,df,11)
      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))stop
     *'root must be bracketed in rtsafe'
      if(fl.eq.0.)then
        rtsafe=x1
        return
      else if(fh.eq.0.)then
        rtsafe=x2
        return
      else if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2)
      dxold=abs(x2-x1)
      dx=dxold
      call funcd(rtsafe,f,df,12)
      do 11 j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0..or. abs(2.*
     *f).gt.abs(dxold*df) ) then
          dxold=dx
          dx=0.5*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return
        call funcd(rtsafe,f,df,13)
        if(f.lt.0.) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
      stop 'rtsafe exceeding maximum iterations'

      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fujl(x,y1,y2,iflag)

c     Evaluate ujl using the hierarchy.


      implicit double precision(a-h,o-z)
      parameter (BIGNO=1.0d10, BIGNI=1.0d-10)

      include 'cmbfast.inc'

c      parameter (m0max=l0max+int(sqrt(dble(40.0*l0max)))+10)
      parameter (m0max=3*l0max+20)
      double precision ak2(m0max),ak2m1(m0max), a2jp1(m0max)
      double precision ujl,ujlp1

      double precision x,y1,y2

      common /param/ beta,ll,K
      common /out/ yout1,yout2

      beta2=beta**2
      beta1=sqrt(beta2-dble(K))
      betam1=1.0d0/beta
      beta1m1=1.0d0/beta1


c Calculating y1,y2 at chi0


c      m=ll+10
      if ((beta.le.ll).and.(K.eq.1)) then
         write(*,*)'Illegal for closed Universes'
         write(*,*)'must have beta> l'
         stop
      end if
      sh=exp(dble(x))/beta
      if (K.eq.-1) then
         m=ll+int(sqrt(dble(40.0*ll)))+10
         chi0=log(sh+sqrt(sh**2+1.0d0))
         ch=sqrt(1.0d0+sh*sh)
      else
         m=ll+int(sqrt(dble(40.0*ll)))+10
         m2=beta-1
         m=min(m,m2)
         chi0=asin(sh)
         ch=sqrt(1.0d0-sh*sh)
      end if
      if (sh.ne.0.0d0)cth=ch/sh

      do lvar=1,m+1
         a2jp1(lvar)=2.0d0*dble(lvar)+1.0d0
      end do


      do lvar=1,m+1
         ak2(lvar)=sqrt(beta2-dble(K*lvar*lvar))
         if (ak2(lvar).gt.0.0d0)ak2m1(lvar)=1.0d0/ak2(lvar)
      end do
      do lvar=m+2,m0max
         ak2(lvar)=0.0d0
      end do

c     setting initial conditions

 10   continue

      cn=cos(beta*chi0)
      sn=sin(beta*chi0)
      uj1=(sn*ch-beta*sh*cn)*betam1*beta1m1/sh**2

      bjp=0.0d0
      bj=1.0d0


      do j3=m,ll+2,-1
         bjm=(a2jp1(j3)*cth*bj-ak2(j3+1)*bjp)*ak2m1(j3)
         bjp=bj
         bj=bjm
         if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
         endif
      end do

      if (((ll+1).gt.m).and.(K.eq.1)) then
         ujlp1=0.0d0
      else
         j3=ll+1
         bjm=(a2jp1(j3)*cth*bj-ak2(j3+1)*bjp)*ak2m1(j3)
         bjp=bj
         bj=bjm
         ujlp1=bjp
      end if

      if ((ll.gt.m).and.(K.eq.1)) then
         ujlp=0.0d0
      else
         j3=ll
         bjm=(a2jp1(j3)*cth*bj-ak2(j3+1)*bjp)*ak2m1(j3)
         bjp=bj
         bj=bjm
         ujl=bjp
      end if

      do j3=ll-1,1,-1
         bjm=(a2jp1(j3)*cth*bj-ak2(j3+1)*bjp)*ak2m1(j3)
         bjp=bj
         bj=bjm
         if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            ujl=ujl*BIGNI
            ujlp1=ujlp1*BIGNI
         endif
      end do

      xnorm=uj1/bjp
      ujl=ujl*xnorm
      ujlp1=ujlp1*xnorm


      bl0=sqrt(beta2-dble(K)*(dble(ll)+1.0d0)**2)


      yout1=ujl
      yout2=(dble(ll)+1.0d0)*ch*ujl-sh*bl0*ujlp1


      y2=dble((yout2-ujl*ch)/yout1/ch)


      argument=yout1*dble(ll)/1.0d-6

c     NEW

      y1=(argument-1.0d0)

c     OLD
cxxx      if (argument.gt.0.0) then
cxxx         y1=log(argument)
cxxx      else
cxxx         write(*,*)'## Error1'
cxxx         write(50,*)beta,ll,chi0
cxxx         write(50,*)ujl,ujlp1
cxxx         write(50,*)x,iflag,ll
cxxx         stop
cxxx      end if


      return
      end
