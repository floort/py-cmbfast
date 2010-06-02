      program jlens
c     Calculates tables for computing Bessel functions in Lensing code.
      implicit double precision (a-h,o-z)

      include 'cmbfast.inc'

      parameter (pi=3.141592654)
      parameter (nx=7*l0max+2)

      double precision x(nx)
      double precision aj0(nx),daj0(nx)
      double precision aj2(nx),daj2(nx)
      double precision aj4(nx),daj4(nx)
      double precision aj6(nx),daj6(nx)
      character*500 filename,arg

      logical intflg

#ifdef IARGC
      if (iargc() .eq. 0) then
         intflg = .true.
      else
         intflg = .false.
         if (iargc() .ne. 2) then
            print *, 'Usage: jlens [lmax] [output filename]'
            stop 1
         end if
         call getarg(1,arg)
         read (arg,*) lmaxim
         call getarg(2,filename)
      endif
#else
      intflg = .true.
#endif

      do i=1,nx
         aj0(i)=0.0d0
         aj2(i)=0.0d0
         aj4(i)=0.0d0
      end do

      if (intflg) then
         write(*,*)'Enter lmax (1500)'
         read(*,*) lmaxim
      endif

c 300 more ls needed for lensing

      lmaxim=lmaxim+300

      xmax=pi*dble(lmaxim)/10.0
      dx=0.05d0
      nt=int(xmax/dx)+2

      if (nt.ge.nx) then
         write(*,*)'You should make nx in jlens.f'
         write(*,*)'and subroutine cllens in lensing.f'
         write(*,*)'bigger than',nt
         write(*,*)'and stard again'
         write(*,*)'The current value is',nx
         stop
      end if

      do i=1,nt
         x(i)=(dble(i)-1.0d0)*dx
         aj0(i)=bessj0(x(i))
         aj2(i)=bessj(2,x(i))
         aj4(i)=bessj(4,x(i))
         aj6(i)=bessj(6,x(i))
      end do

      ajlo=1.0d40
      ajhi=1.0d40
      call spline(x,aj0,nt,ajlo,ajhi,daj0)
      call spline(x,aj2,nt,ajlo,ajhi,daj2)
      call spline(x,aj4,nt,ajlo,ajhi,daj4)
      call spline(x,aj6,nt,ajlo,ajhi,daj6)

      if (intflg) then
         write(*,*)'Enter output file (jlens.dat)'
         read(*,'(a500)')filename
      endif
      open(unit=9,file=filename
     <     ,status='unknown',form='unformatted')
      rewind 9
      write(9)nx
      write(9)lmaxim
      write(9)nt
      write(9)dx
      write(9)aj0
      write(9)daj0
      write(9)aj2
      write(9)daj2
      write(9)aj4
      write(9)daj4
      write(9)aj6
      write(9)daj6
      close(9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION bessj1(x)
      DOUBLE PRECISION  bessj1,x
      DOUBLE PRECISION  ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.0d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=sqrt(.636619772d0/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,real(x))
      endif
      return
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION bessj(n,x)
      INTEGER n,IACC
      DOUBLE PRECISION  bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
CU    USES bessj0,bessj1
      INTEGER j,jsum,m
      DOUBLE PRECISION  ax,bj,bjm,bjp,sum,tox,bessj0,bessj1

      if(n.lt.2)stop 'bad argument n in bessj'
      ax=abs(x)
      if(ax.eq.0.)then
        bessj=0.
      else if(ax.gt.float(n))then
        tox=2.0d0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2.0d0/ax
        m=2*((n+int(sqrt(float(IACC*n))))/2)
        bessj=0.0d0
        jsum=0
        sum=0.0d0
        bjp=0.0d0
        bj=1.0d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.0d0*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0.0d0.and.mod(n,2).eq.1)bessj=-bessj
      return
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
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
