cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C This program calculates and stores the spherical Bessel functions
C for the cmb integrator.
C Developed by Uros Seljak (useljak@cfa.harvard.edu)
C and Matias Zaldarriaga (matiasz@arcturus.mit.edu)
C based on a code written by Arthur Kosowsky   akosowsky@cfa.harvard.edu
C
C Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
C the Massachusetts Institute of Technology.  All rights reserved.
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program jlgen
      include 'cmbfast.inc'

      parameter (xlimmin=10.0)
      real x,ajl
      integer l(lmax),i,j,lmo
      character*500 filename,arg
      logical intflg
      common /lvalues1/ l,l0,lmo
      save /lvalues1/

#ifdef IARGC
      if (iargc() .eq. 0) then
         intflg = .true.
      else
         intflg = .false.
         if (iargc() .ne. 3) then
            print *, 'Usage: jlgen [lmax] [kmax] [output filename]'
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
         write(*,*)'Maximum value of l, keta: (1500, 3000)'
         read(*,*)lmoin,kmax0
      endif
      kmax=kmax0-25+151

c For lensing the spectrum we need some extra ls. So even
c if the user askes for lmo we will calculate to lmo+300

       if ((lmoin+300).gt.l0max) then
          write(*,*)'Sorry there is a maximum l of ', l0max
          write(*,*)'The code needs to calculate up to',lmoin+300
          write(*,*)'300 more than what you asked.'
          write(*,*)'Please increase l0max in both'
          write(*,*)'jlgen.f, ujlgen.f,'
          write(*,*)'jlens.f,ccmbflat.f and cmbopen.f'
          write(*,*)'in all places where it appears'
          stop
       end if
       call initlval(lmoin)

        if (intflg) then
           write(*,*) 'Enter output filename'
           read(*,'(a500)') filename
        endif
        open(unit=9,file=filename
     <,status='unknown',form='unformatted')
        rewind 9
        write(9)lmo
        write(9)kmax0

c     For safety purposes check. We will use this to
c     check if the lvalues are the same in file and code.
        write(9)l0
        do j=1,l0
           write(9)l(j)
        end do

       do 30 j=1,l0
        do 40 i=1,kmax
           if (i.le.151) then
           if (i.le.51) then
              x=(i-1)/10.0
              else
              x=(i-51)/5.0+5.0
           end if
           else
              x=(i-151)+25.0
           end if
         xlim=l(j)/20.0
         xlim=max(xlim,xlimmin)
         xlim=l(j)-xlim
         if (x-xlim.gt.1E-10) then
          call bjl(l(j),x,ajl)
cxxx          if ((l(j).eq.3).and.(x.le.0.2)) then
cxxx             ajl=0.0
cxxx          end if
cxxx          if ((l(j).gt.3).and.(x.lt.0.6)) then
cxxx             ajl=0.0
cxxx          end if
cxxx          if ((l(j).ge.5).and.(x.lt.1.0)) then
cxxx             ajl=0.0
cxxx          end if
          write(9)ajl
         endif
 40     continue
30     continue
       stop
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Subroutine version of the main jlgen program.
        subroutine(lmoin,kmax0,filename)
      include 'cmbfast.inc'

      parameter (xlimmin=10.0)
      real x,ajl
      integer l(lmax),i,j,lmo
      character*500 filename,arg
      common /lvalues1/ l,l0,lmo
      save /lvalues1/

      kmax=kmax0-25+151

c For lensing the spectrum we need some extra ls. So even
c if the user askes for lmo we will calculate to lmo+300

       if ((lmoin+300).gt.l0max) then
          stop
       end if
       call initlval(lmoin)

        open(unit=9,file=filename
     <,status='unknown',form='unformatted')
        rewind 9
        write(9)lmo
        write(9)kmax0

c     For safety purposes check. We will use this to
c     check if the lvalues are the same in file and code.
        write(9)l0
        do j=1,l0
           write(9)l(j)
        end do

       do 30 j=1,l0
        do 40 i=1,kmax
           if (i.le.151) then
           if (i.le.51) then
              x=(i-1)/10.0
              else
              x=(i-51)/5.0+5.0
           end if
           else
              x=(i-151)+25.0
           end if
         xlim=l(j)/20.0
         xlim=max(xlim,xlimmin)
         xlim=l(j)-xlim
         if (x-xlim.gt.1E-10) then
          call bjl(l(j),x,ajl)
cxxx          if ((l(j).eq.3).and.(x.le.0.2)) then
cxxx             ajl=0.0
cxxx          end if
cxxx          if ((l(j).gt.3).and.(x.lt.0.6)) then
cxxx             ajl=0.0
cxxx          end if
cxxx          if ((l(j).ge.5).and.(x.lt.1.0)) then
cxxx             ajl=0.0
cxxx          end if
          write(9)ajl
         endif
 40     continue
30     continue
       stop
       end




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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine bjl(L,x,jl)
c  Calculates the spherical bessel function j_l(x)

c  and optionally its derivative for real x and integer l>=0.

c  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from

c  G.N.Watson, A Treatise on the Theory of Bessel Functions,
c  2nd Edition (Cambridge University Press, 1944).
c  Higher terms in expansion for x near l given by
c  Airey in Phil. Mag. 31, 520 (1916).

c  This approximation is accurate to near 0.1% at the boundaries
c  between the asymptotic regions; well away from the boundaries
c  the accuracy is better than 10^{-5}. The derivative accuracy
c  is somewhat worse than the function accuracy but still better
c  than 1%.

c  Point *jlp initially to a negative value to forego calculating
c  the derivative; point it to a positive value to do the derivative
c  also (Note: give it a definite value before the calculation
c  so it's not pointing at junk.) The derivative calculation requires

c  only arithmetic operations, plus evaluation of one sin() for the
c  x>>l region.


c  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu
c  This fortran version only computes j_l(x)
        implicit double precision(a-h,o-z)
        integer L
        double precision nu, nu2,ax,ax2,beta,beta2,beta4,beta6
        double precision sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1
        double precision cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b
        double precision trigarg,trigcos,expterm,prefactor,llimit,
     &                   ulimit,fl
        real x,jl

        PI = 3.1415926536
        ROOTPI = 1.772453851
        GAMMA1 = 2.6789385347             !/* Gamma function of 1/3 */
        GAMMA2 = 1.3541179394             !/* Gamma function of 2/3 */


        ax = abs(dble(x))
        fl = l


        beta = fl**0.325
        llimit=1.31*beta   !/* limits of asymptotic regions; fitted */
        ulimit=1.48*beta

         nu= fl + 0.5

         nu2=nu*nu

        if (l .lt. 0) then
                print*, 'Bessel function index < 0\n'
                stop
        end if

c          /************* Use closed form for l<6 **********/

        if (l .lt. 6) then

        sx=sin(ax)
        cx=cos(ax)
        ax2=ax*ax

            if(l .eq. 0) then
                 if(ax .gt. 0.001) then
                 jl=real(sx/ax)
                 else

                 jl=real(1.0d0-ax2/6.0d0)
                 end if                         !   /* small x trap */
            endif


            if(l .eq. 1) then

                 if(ax .gt. 0.001) then
                 jl=real((sx/ax -cx)/ax)
                 else

                 jl=real(ax/3.0d0)
                 end if
            endif

            if(l .eq. 2) then
                  if(ax .gt. 0.001) then
                  jl=real((-3.0d0*cx/ax
     1                    -sx*(1.0d0-3.0d0/ax2))/ax)
                  else

                  jl=real(ax2/15.0d0)
                  end if
            endif

            if(l .eq. 3) then
                if(ax .gt. 0.001) then
                jl=real((cx*(1.0d0-15.0d0/ax2)
     1                    -sx*(6.0d0-15.0d0/ax2)/ax)/ax)
                else

                jl=real(ax*ax2/105.0d0)
                endif
            endif

            if(l .eq. 4) then
                if(ax .gt. 0.001) then
        jl=real((sx*(1.0d0-45.0d0/(ax*ax)+105.0d0
     1       /(ax*ax*ax*ax)) +cx*(10.0d0-105.0d0/(ax*ax))/ax)/ax)
                else

                jl=real(ax2*ax2/945.0d0)
                end if
            endif

             if(l .eq. 5) then

                if(ax .gt. 0.001) then
        jl=real((sx*(15.0d0-420.0d0/(ax*ax)+945.0d0
     1    /(ax*ax*ax*ax))/ax -cx*(1.0-105.0d0/(ax*ax)+945.0d0
     1                                    /(ax*ax*ax*ax)))/ax)
                else

                jl=real(ax2*ax2*ax/10395.0d0)
                endif
             endif


c          /********************** x=0 **********************/

        else if (ax .lt. 1.d-30) then
        jl=0.0

c          /*************** Region 1: x << l ****************/

        else if (ax .le. fl+0.5-llimit) then


c       beta=acosh(nu/ax)
        if (nu/ax .lt. 1.d0) print*, 'trouble with acosh'
        beta = dlog(nu/ax + dsqrt((nu/ax)**2 - 1.d0) )
                !(4.6.21)
        cotb=nu/sqrt(nu*nu-ax*ax)      ! /* cotb=coth(beta) */
        cot3b=cotb*cotb*cotb
        cot6b=cot3b*cot3b
        secb=ax/nu
        sec2b=secb*secb
        sec4b=sec2b*sec2b
        sec6b=sec4b*sec2b
        sum1=2.0+3.0*sec2b
        expterm=sum1*cot3b/(24.0*nu)
        sum2=4.0+sec2b
        expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2)
        sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
        expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2)
        sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
        expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
        expterm=exp(-nu*beta+nu/cotb-expterm)
        prefactor=sqrt(cotb/secb)/(2.0*nu)
        jl=real(prefactor*expterm)

c          /**************** Region 2: x >> l ****************/


        else if (ax .ge. fl+0.5+ulimit) then


        beta=acos(nu/ax)
        cotb=nu/sqrt(ax*ax-nu*nu)      !/* cotb=cot(beta) */
        cot3b=cotb*cotb*cotb
        cot6b=cot3b*cot3b
        secb=ax/nu
        sec2b=secb*secb
        sec4b=sec2b*sec2b
        sec6b=sec4b*sec2b
        trigarg=nu/cotb - nu*beta - PI/4.0
        sum1=2.0+3.0*sec2b
        trigarg = trigarg - sum1*cot3b/(24.0*nu)
        sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b
        trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2)
        trigcos=cos(trigarg)
        sum2=4.0+sec2b
        expterm=sum2*sec2b*cot6b/(16.0*nu2)
        sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b
        expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2)
        expterm=exp(-expterm)
        prefactor=sqrt(cotb/secb)/nu
        jl=real(prefactor*expterm*trigcos)

c          /***************** Region 3: x near l ****************/

        else



        beta=ax-nu

        beta2=beta*beta
        beta4=beta2*beta2
        beta6=beta2*beta4
        sx=6.0/ax
        sx2=sx*sx
        cx=sqrt(sx)

        secb=sx**0.333333333

        sec2b=secb*secb

        deriv1=GAMMA1*secb
        deriv1= deriv1+ beta*GAMMA2*sec2b
        sum1=(beta2/6.0-1.0/15.0)*beta
        deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0
        sum2=beta4/24.0-beta2/24.0+1.0/280.0
        deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0
        sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0
        deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0
        sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0-
     2  13.0/31500.0)*beta
        deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0
        sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0
     1               -121.0*beta2/907200.0 + 7939.0/232848000.0)*beta
        deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0

        jl=real(deriv1*cx/(12.0*ROOTPI))

        end if

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
