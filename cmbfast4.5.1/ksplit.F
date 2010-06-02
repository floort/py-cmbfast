C Integrator for CMB anisotropy, CMB polarization and transfer functions
C Developed by Uros Seljak (useljak@princeton.edu)
C and Matias Zaldarriaga (matiasz@sns.ias.edu).
C See the LICENSE file for restrictions on the use, modification and
C distribution of this software.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This code is a wrapper to call cmbfast using the k-split method
c and to combine the resulting spectra.
c
        program ksplit
        implicit double precision(a-h,o-z)

        integer slopes, type, lmax

        real parvec(7), nvec(20),alvec(20), ketamax
        character*500 outfile

        double precision wdyn
        integer ndyn
        common /qparm/ wdyn,ndyn

        write(*,*)' '
        write(*,*)' '
        write(*,*)'This is a driver that illustrates',
     $       'k-splitting method.'
        write(*,*)'Bessel functions need to be precomputed'
        write(*,*)'and stored in jl.dat and ujl.dat'
        write(*,*)'should be used with care and accuracy should'
        write(*,*)'be checked. The k-split method is particularly'
        write(*,*)'useful when running grids. In that case you'
        write(*,*)'will need to create your own driver'
        write(*,*)' '
        write(*,*)' '

        write (*,*) 'Enter(0) for no quintessence models '
        write(*,*)'(1) constant equation of state'
        read(*,*)ndyn


c to run other models you have to modify subroutine gentable acordingly

        if(ndyn.eq.1)then
          write (*,*) 'Enter wdyn (e.g -1)'
          read (*,*) wdyn
        else
           wdyn=-1.0
        endif

        write(*,*)'Enter lmax,ketamax for the high model'
        read(*,*)lmax,ketamax


        write(*,*)'Enter tau, Omega_k, Omega_de,
     1    omega_dm, omega_b, f_nu, ns'
        read(*,*)parvec

        do i=1,20
           nvec(i)=0.0d0
           alvec(i)=0.0d0
        end do
        nvec(1)=parvec(7)
        slopes=1
        tau=parvec(1)
        if ((parvec(1).eq.0.0).and.
     $       (parvec(2).eq.0.0).and.(parvec(3).eq.0.0)) then
           write(*,*)'The model is flat with no Lambda or tau'
           write(*,*)'will not do the k-split'
           outfile='cl.dat'

           call computeCl(parvec,1,nvec,alvec,0,lmax,ketamax,outfile)

        else

c       Compute low ls
           outfile='cl_low.dat'

           write(*,*)'1. output file ..... ',outfile

           call computeCl(parvec,1,nvec,alvec,2,400,400.0,outfile)


c       Compute high ls
           parvec(1)=0.0
           parvec(2)=0.0
           parvec(3)=0.0
           outfile='cl_high.dat'

           call computeCl(parvec,1,nvec,alvec,3,lmax,ketamax,outfile)

c       Combine spectra
           call combineCl(tau)
        end if

        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine computeCl(parvec,slopes,nvec,alvec,
     $     type,lmax,ketamax,outfile)

        implicit none
        real Omega_b,Omega_c,Omega_v,Omega_n,Omega_m,Omega_d,f_nu
        real H_0,t_cmb,Y_he,a_nu_r,a_nu_nr,tau
        real Omega_k, h,rwdyn

        character*500 outfile
        real parvec(7), nvec(20), alvec(20),ketamax
        integer slopes, type, lmax, i
        double precision wdyn

        integer ndyn

!       Parameter       Meaning
!       1               tau
!       2               Omega_k
!       3               Omega_Lambda
!       4               Omega_dm h^2 = (Omega_cdm + Omega_nu) h^2
!       5               Omega_baryon h^2
!       6               f_nu = Omega_nu/Omega_dm
!       7               ns or nt

        tau     = parvec(1)
        Omega_k = parvec(2)

        Omega_v = parvec(3)
        omega_m = parvec(4) + parvec(5)
        f_nu    = parvec(6)
#ifdef DIM
        omegamatter=(1.0-parvec(2))-2*(1-parvec(2))) -parvec(3)
        if (1.-Omega_k-Omega_v.le.0.) pause 'DEATH ERROR: Omega_m NEGATIVE'
        h       = sqrt(omega_m/omegamatter)
#else

        if (1.-Omega_k-Omega_v.le.0.)
     $stop 'DEATH ERROR: Omega_m NEGATIVE'
        h       = sqrt(omega_m/(1.-Omega_k-Omega_v))
#endif
        Omega_d = parvec(4)/h**2                ! (h^2 Od)/h^2
        omega_b = parvec(5)/h**2                ! (h^2 Ob)/h^2
        Omega_c = Omega_d*(1.-f_nu)
        Omega_n = Omega_d*f_nu
        H_0     = 100.*h
        t_cmb   = 2.726
        Y_he    = 0.24
        if (abs(Omega_n).lt.1e-5) then
          a_nu_r        = 3.
          a_nu_nr       = 0.
        else
          a_nu_r        = 2.
          a_nu_nr       = 1.
        end if

        rwdyn=-1.0

c   TESTING
        print *,'CMBFAST called with...'
        write(*,'(1f6.1,1i5,2i4,99f8.3)') ketamax, lmax, type,
     &       slopes, tau,
     &       Omega_k, Omega_v, Omega_c, Omega_b, Omega_n, h,rwdyn

        call run_cmbfast(
     &          Omega_b,Omega_c,Omega_v,Omega_n,rwdyn,
     &          H_0,t_cmb,Y_he,a_nu_r,a_nu_nr,
     &          tau,slopes,nvec,alvec,type,lmax,ketamax,outfile)
        return
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine combineCl(tau)

        implicit double precision(a-h,o-z)
        include 'cmbfast.inc'

        character*500 temp

        double precision llow(l0max),ctlow(l0max),celow(l0max),
     &                   cclow(l0max)
        double precision lhig(l0max),cthig(l0max),cehig(l0max),
     &                   cchig(l0max)
        double precision cthigpr(l0max),cehigpr(l0max),cchigpr(l0max)

c       Read low model

        open(unit=3,file='cl_low.dat',form='formatted'
     <,status='old')
        read(3,*)temp
        read(3,*)temp
        read(3,*)dlsslow
        i=0
 5      continue
        read(3,*,end=10)l,ct,ce,cc
        i=i+1
        llow(i)=dble(l)
        ctlow(i)=ct
        celow(i)=ce
        cclow(i)=cc
        goto 5
 10     continue
        nlow=i
        lmaxlow=int(llow(nlow))
        close(3)


c       Read high model

        open(unit=3,file='cl_high.dat',form='formatted'
     <,status='old')
        read(3,*)temp
        read(3,*)temp
        read(3,*)dlsshig
        i=0
 15     continue
        read(3,*,end=20)l,ct,ce,cc
        i=i+1
        lhig(i)=dble(l)
        cthig(i)=ct
        cehig(i)=ce
        cchig(i)=cc
        goto 15
 20     continue
        close(3)
        nhigh=i

c       Shift in l and tau correction
        do i=1,nhigh
           lhig(i)=lhig(i)*dlsslow/dlsshig
           cthig(i)=cthig(i)*exp(-2.0*tau)
           cehig(i)=cehig(i)*exp(-2.0*tau)
           cchig(i)=cchig(i)*exp(-2.0*tau)
        end do
        lmaxhig=int(lhig(nhigh))
        lminhig=int(lhig(1))+1


c       Spline high model back to uniform l.
        call spline(lhig,cthig,nhigh,1.0d40,1.0d40,cthigpr)
        call spline(lhig,cehig,nhigh,1.0d40,1.0d40,cehigpr)
        call spline(lhig,cchig,nhigh,1.0d40,1.0d40,cchigpr)

        open(unit=3,file='cl_combine.dat',form='formatted'
     <,status='unknown')

        llo=1

        do il=2,lmaxhig
           ctint=0.0d0
           ceint=0.0d0
           ccint=0.0d0
           xi=dble(il)
           if (il.gt.lminhig) then

              templlo=xi/dlsslow*dlsshig-1.0d0
              llo=int(templlo)
              lhi=llo+1
              ho=lhig(lhi)-lhig(llo)
              a0=(lhig(lhi)-xi)/ho
              b0=(xi-lhig(llo))/ho
              ctint=a0*cthig(llo)+b0*cthig(lhi)+((a0**3-a0)*
     1             cthigpr(llo)+(b0**3-b0)
     2             *cthigpr(lhi))*ho*ho /6.d0
              ceint=a0*cehig(llo)+b0*cehig(lhi)+((a0**3-a0)
     1             *cehigpr(llo)+(b0**3-b0)
     2             *cehigpr(lhi))*ho*ho /6.d0
              ccint=a0*cchig(llo)+b0*cchig(lhi)+((a0**3-a0)
     1             *cchigpr(llo)+(b0**3-b0)
     2             *cchigpr(lhi))*ho*ho /6.d0
           end if

           if (il.lt.lmaxlow)then
              ctint=ctint+ctlow(il-1)
              ceint=ceint+celow(il-1)
              ccint=ccint+cclow(il-1)
           end if

c       TESTING
c          write(16,'(1I6,200E15.5)')llo,xi,lhig(llo),
c       1       lhig(lhi),a0,b0,ho,ctint
c
c          write(17,'(200E15.5)')xi,cthig(llo),cthig(lhi),
c       1       cthigpr(llo),cthigpr(lhi)

           write(3,'(1I5,3E13.5)')il,ctint,ceint,ccint

        end do
        close(3)

        return
        end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
