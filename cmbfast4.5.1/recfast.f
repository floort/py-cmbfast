cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Integrator for Cosmic Recombination of Hydrogen and Helium,
C developed by Douglas Scott (dscott@astro.ubc.ca)
C based on calculations in the papers Seager, Sasselov & Scott
C (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000).
C
C Permission to use, copy, modify and distribute without fee or royalty at
C any tier, this software and its documentation, for any purpose and without
C fee or royalty is hereby granted, provided that you agree to comply with
C the following copyright notice and statements, including the disclaimer,
C and that the same appear on ALL copies of the software and documentation,
C including modifications that you make for internal use or for distribution:
C
C Copyright 1999 by University of British Columbia.  All rights reserved.
C
C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
C REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  
C BY WAY OF EXAMPLE, BUT NOT LIMITATION,
c U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF 
C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE 
C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CN	Name:	RECFAST
CV	Version: 1.2
C 
CP	Purpose:  Calculate ionised fraction as a function of redshift.
CP		  Solves for H and He simultaneously, and includes
CP		  "fudge factor" for low z effect.
C
CD	Description: Solves for ionisation history since recombination
CD	using the equations in Seager, Sasselov & Scott (ApJ, 1999).
CD	The Cosmological model can be flat or open.
CD	The matter temperature is also followed.
CD	The values for \alpha_B for H are from Hummer (1994).
CD	The singlet HeI coefficient is a fit from the full code.
CD	Care is taken to use the most accurate constants.
CD	Note that some parameters are fixed (e.g. N_nu=3, nu's are
CD	massless, w=-1, etc.) - some users may want to explictly
CD	imput their own H(z) to account for extra physics.
CD	This is provided as a PROGRAM, which can be easily converted
CD	to a SUBROUTINE for use in CMB Boltzmann codes.
C		
CA	Arguments:
CA	Name, Description
CA	Double precision throughout
CA
CA	z is redshift - W is sqrt(1+z), like conformal time
CA	x is total ionised fraction, relative to H
CA	x_H is ionized fraction of H - y(1) in R-K routine
CA	x_He is ionized fraction of He - y(2) in R-K routine
CA	Tmat is matter temperature - y(3) in R-K routine
CA	f's are the derivatives of the Y's
CA	alphaB is case B recombination rate
CA	alpHe is the singlet only HeII recombination rate
CA	a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	a_VF is Verner and Ferland type fitting parameter for Helium
CA	b_VF is Verner and Ferland type fitting parameter for Helium
CA	T_0 is Verner and Ferland type fitting parameter for Helium
CA	T_1 is Verner and Ferland type fitting parameter for Helium
CA	Tnow is the observed CMB temperature today
CA	OmegaT is the total Omega_0
CA      OmegaL is the Omega_0 contribution from a Cosmological constant
CA      OmegaK is the Omega_0 contribution in curvature (1-O_T-O_L)
CA      OmegaB is Omega in baryons today
CA	OmegaC is the Omega_0 in (cold) dark matter: OmegaT=OmegaC+OmegaB
CA	Yp is the primordial helium abundace
CA	fHe is He/H number ratio = Yp/4(1-Yp)
CA	Trad and Tmat are radiation and matter temperatures
CA	OmegaB is Omega in baryons today
CA	H is Hubble constant in units of 100 km/s/Mpc
CA	HOinp is input value of Hubble constant in units of 100 km/s/Mpc
CA	HO is Hubble constant in SI units
CA	bigH is 100 km/s/Mpc in SI units
CA	G is grvitational constant
CA	n is number density of hydrogen
CA	Nnow is number density today
CA	x0 is initial ionized fraction
CA	x_H0 is initial ionized fraction of Hydrogen
CA	x_He0 is initial ionized fraction of Helium
CA	rhs is dummy for calculating x0
CA	zinitial and zfinal are starting and ending redshifts
CA	fnu is the contribution of neutrinos to the radn. energy density
CA	zeq is the redshift of matter-radiation equality
CA	zstart and zend are for each pass to the integrator
CA	w0 and w1 are conformal-time-like initial and final zi and zf's
CA	Lw0 and Lw1 are logs of w0 and w1
CA	hw is the interval in W
CA	C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
CA	m_e,m_H: electron mass and mass of H atom in SI
CA	sigma: Thomson cross-section
CA	a: radiation constant for u=aT^4
CA	Pi: Pi
CA	Lambda: 2s-1s two photon rate for Hydrogen
CA	Lambda_He: 2s-1s two photon rate for Helium
CA	DeltaB: energy of first excited state from continuum = 3.4eV
CA	DeltaB_He: energy of first excited state from cont. for He = 3.4eV
CA	L_H_ion: level for H ionization in m^-1
CA	L_H_alpha: level for H Ly alpha in m^-1
CA	L_He1_ion: level for HeI ionization
CA	L_He2_ion: level for HeII ionization
CA	L_He_2s: level for HeI 2s
CA	L_He_2p: level for HeI 2p
CA	Lalpha: Ly alpha wavelength in SI
CA	Lalpha_He: Helium I 2p-1s wavelength in SI
CA	mu_H,mu_T: mass per H atom and mass per particle
CA	H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
CA	CDB=DeltaB/k_B			Constants derived from B1,B2,R
CA	CDB_He=DeltaB_He/k_B		n=2-infinity for He in Kelvin
CA	CB1=CDB*4.			Lalpha and sigma_Th, calculated
CA	CB1_He1: CB1 for HeI ionization potential
CA	CB1_He2: CB1 for HeII ionization potential
CA	CR=2*Pi*(m_e/h_P)*(k_B/h_P)	once and passed in a common block
CA	CK=Lalpha**3/(8.*Pi)
CA	CK_He=Lalpha_He**3/(8.*Pi)
CA	CL=C*h_P/(k_B*Lalpha)
CA	CL_He=C*h_P/(k_B*Lalpha_He)
CA	CT=(8./3.)*(sigma/(m_e*C))*a
CA	Bfact=exp((E_2p-E_2s)/kT)	Extra Boltzmann factor
CA	tol: tolerance for the integrator
CA	cw(24),w(3,9): work space for DVERK
CA	Ndim: number of d.e.'s to solve (integer)
CA	Nz: number of output redshitf (integer)
CA	I: loop index (integer)
CA	ind,nw: work-space for DVERK (integer)
C
CG	Global data (common blocks) referenced:
CG	/zLIST/zinitial,zfinal,Nz
CG	/Cfund/C,k_B,h_P,m_e,m_H,sigma,a,Pi
CG	/data/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
CG		fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He
CG      /Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C
CF	File & device access:
CF	Unit	/I,IO,O	/Name (if known)
C
CM	Modules called:
CM	DVERK (numerical integrator)
CM	GET_INIT (initial values for ionization fractions)
CM	ION (ionization and Temp derivatices)
C
CC	Comments:
CC	none
C
CH	History:
CH	CREATED		(simplest version) 19th March 1989
CH	RECREATED	11th January 1995
CH			includes variable Cosmology
CH			uses DVERK integrator
CH			initial conditions are Saha
CH	TESTED		a bunch, well, OK, not really
CH	MODIFIED	January 1995 (include Hummer's 1994 alpha table)
CH			January 1995 (include new value for 2s-1s rate)
CH			January 1995 (expand comments)
CH			March 1995 (add Saha for Helium)
CH			August 1997 (add HeII alpha table)
CH			July 1998 (include OmegaT correction and H fudge factor)
CH			Nov 1998 (change Trad to Tmat in Rup)
CH			Jan 1999 (tidied up for public consumption)
CH			Sept 1999 (switch to formula for alpha's, fix glitch)
CH			Feb 2000 (fixed overflow problem in He_Boltz)
CH			Oct 2001 (fixed OmegaT in z_eq)
CH			June 2003 (fixed error in Rdown_He formula)
CH			June 2003 (fixed He recombination coefficients)
CH			June 2003 (comments to point out fixed N_nu etc.)
CH				
C-
C	===============================================================

        subroutine recfast(OmegaB,OmegaC,omegav,HOinp,tcmb,Yp,annur)

	implicit none

	integer Ndim,Nz,I,Nz0
        parameter (Nz0=10000)
C	--- Arguments
        double precision zrec(Nz0),xrec(Nz0),dxrec(Nz0)
	double precision Trad,Tmat,d0hi,d0lo
        double precision OmegaT,OmegaB,H,HO,HOinp,bigH,G,
     $       OmegaL,OmegaK,OmegaC
	double precision z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0,annur,tcmb
	double precision Tnow,zinitial,zfinal,Nnow,z_eq,fnu,omegav
	double precision zstart,zend,w0,w1,Lw0,Lw1,hw
	double precision C,k_B,h_P,m_e,m_H,sigma,a,Pi
	double precision Lambda,DeltaB,DeltaB_He,Lalpha,
     $       mu_H,mu_T,H_frac
	double precision Lambda_He,Lalpha_He,Bfact,CK_He,CL_He
	double precision L_H_ion,L_H_alpha,L_He1_ion,L_He2_
     $       ion,L_He_2s,L_He_2p
	double precision CB1,CDB,CR,CK,CL,CT,Yp,fHe,CB1_He1,
     $       CB1_He2,CDB_He,fu
	double precision tol
	double precision cw(24),w(3,9)
	double precision y(3)

	integer ind,nw

	character*80 fileout

C	--- Parameter statements
	parameter(bigH=100.0D3/(1.0D6*3.0856775807D16))	!Ho in s-1
	parameter(tol=1.D-5)				!Tolerance for R-K

	external ION

C	--- Commons
c	common/zLIST/zinitial,zfinal,Nz
	common/Cfund/C,k_B,h_P,m_e,m_H,sigma,a,Pi
	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,
     $       CK_He,CL_He,fu
	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
        common/recout/zrec,xrec,dxrec,Nz
        save/recout/

C	===============================================================

C	--- Data
	data	C,k_B,h_P	/2.99792458D8,1.380658D-23,6.6260755D-34/
	data	m_e,m_H		/9.1093897D-31,1.673725D-27/	!av. H atom
	data	sigma,a		/6.6524616D-29,7.565914D-16/
	data	Pi		/3.141592653589d0/
	data	G		/6.67259D-11/
C	Fundamental constants in SI units

	data	Lambda		/8.2245809d0/
	data	Lambda_He	/51.3d0/	!new value from Dalgarno
	data	L_H_ion		/1.096787737D7/	!level for H ion. (in m^-1)
	data	L_H_alpha	/8.225916453D6/ !averaged over 2 levels
	data	L_He1_ion	/1.98310772D7/	!from Drake (1993)
	data	L_He2_ion	/4.389088863D7/	!from JPhysChemRefData (1987)
	data	L_He_2s		/1.66277434D7/	!from Drake (1993)
	data	L_He_2p		/1.71134891D7/	!from Drake (1993)
C	2 photon rates and atomic levels in SI units

C	dimensions for integrator
	Ndim = 3

c	write(*,*)'recfast version 1.2'
c	write(*,*)'Using Hummer''s case B recombination rates for H'
c	write(*,*)' with fudge factor = 1.14,'
c	write(*,*)' and a fit to tabulated HeII singlet recombination rates'
c	write(*,*)

c	These are easy to inquire as input, but let's use simple values
	zinitial = 1.d4
	z = zinitial
	zfinal=0.d0
c	will output every 10 in z, but this is easily changed also

	OmegaT=OmegaC+OmegaB            !total dark matter + baryons
	OmegaK=1.d0-OmegaT-OmegaL	!curvature
        Omegal=omegav
        Tnow=tcmb

c	convert the Hubble constant units
	H = HOinp/100.d0
	HO = H*bigH


c	sort out the helium abundance parameters
	mu_H = 1.d0/(1.d0-Yp)			!Mass per H atom
	mu_T = 4.d0/(4.d0-3.d0*Yp)		!Mass per atom
	fHe = Yp/(4.d0*(1.d0-Yp))		!n_He_tot / n_H_tot

	Nnow = 3.d0*HO*HO*OmegaB/(8.d0*Pi*G*mu_H*m_H)
	n = Nnow * (1.d0+z)**3
	fnu = (annur*7.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
	z_eq = (3.d0*(HO*C)**2/(8.d0*Pi*G*a*(1.d0+fnu)*Tnow**4))*OmegaT
	z_eq = z_eq - 1.d0

C	Set up some constants so they don't have to be calculated later
	Lalpha = 1.d0/L_H_alpha
	Lalpha_He = 1.d0/L_He_2p
	DeltaB = h_P*C*(L_H_ion-L_H_alpha)
	CDB = DeltaB/k_B
	DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)	!2s, not 2p
	CDB_He = DeltaB_He/k_B
	CB1 = h_P*C*L_H_ion/k_B
	CB1_He1 = h_P*C*L_He1_ion/k_B	!ionization for HeI
	CB1_He2 = h_P*C*L_He2_ion/k_B	!ionization for HeII
	CR = 2.d0*Pi*(m_e/h_P)*(k_B/h_P)
	CK = Lalpha**3/(8.d0*Pi)
	CK_He = Lalpha_He**3/(8.d0*Pi)
	CL = C*h_P/(k_B*Lalpha)
	CL_He = C*h_P/(k_B/L_He_2s)	!comes from det.bal. of 2s-1s
	CT = (8.d0/3.d0)*(sigma/(m_e*C))*a
	Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

C	Matter departs from radiation when t(Th) > H_frac * t(H)
C	choose some safely small number
	H_frac = 1.D-3

c	Fudge factor to approximate for low z out of equilibrium effect
	fu=1.14d0

c	Set initial matter temperature
	y(3) = Tnow*(1.d0+z)            !Initial rad. & mat. temperature
	Tmat = y(3)

	call get_init(z,x_H0,x_He0,x0)

	y(1) = x_H0
	y(2) = x_He0

c	OK that's the initial conditions, now start writing output file

	w0=1.d0/ dsqrt(1.d0 + zinitial)	!like a conformal time
	w1=1.d0/ dsqrt(1.d0 + zfinal)
	Lw0 = dLog(w0)
	Lw1 = dLog(w1)
	Nz=10000
	hW=(Lw1-Lw0)/dfloat(Nz)		!interval in log of conf time

c	Set up work-space stuff for DVERK
	ind  = 1
	nw   = 3
	do i = 1,24
	  cw(i) = 0.d0
	end do

	do i = 1,Nz
C       calculate the start and end redshift for the interval at each z
C	or just at each z
	  zstart = zinitial + dfloat(i-1)*(zfinal-zinitial)/dfloat(Nz)
	  zend   = zinitial + dfloat(i)*(zfinal-zinitial)/dfloat(Nz)

C Use Saha to get x_e, using the equation for x_e for ionized helium
C and for neutral helium.
C Everything ionized above z=8000.  First ionization over by z=5000.
C Assume He all singly ionized down to z=3500, then use He Saha until
C He is 99% singly ionized, and *then* switch to joint H/He recombination.

	  z = zend

	  if (zend.gt.8000.d0) then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    x0 = 1.d0+2.d0*fHe
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if(z.gt.5000.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He2/(Tnow*(1.d0+z)) ) / Nnow
	    rhs = rhs*1.d0		!ratio of g's is 1 for He++ <-> He+
	    x0 = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2
	1	+ 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if(z.gt.3500.d0)then

	    x_H0 = 1.d0
            x_He0 = 1.d0
            x0 = x_H0 + fHe*x_He0
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

          else if(y(2).gt.0.99)then

	    x_H0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He1/(Tnow*(1.d0+z)) ) / Nnow
	    rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs )
	1	- (rhs-1.d0))
	    x0 = x_He0
	    x_He0 = (x0 - 1.d0)/fHe
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if (y(1).gt.0.99d0) then

	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1/(Tnow*(1.d0+z)) ) / Nnow
	    x_H0 = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )

	    call DVERK(nw,ION,zstart,y,zend,tol,ind,cw,nw,w)
	    y(1) = x_H0
	    x0 = y(1) + fHe*y(2)

	  else

	    call DVERK(nw,ION,zstart,y,zend,tol,ind,cw,nw,w)

	    x0 = y(1) + fHe*y(2)

	  end if

	  Trad = Tnow * (1.d0+zend)
	  Tmat = y(3)
	  x_H = y(1)
	  x_He = y(2)
	  x = x0

          zrec(i)=zend
          xrec(i)=x

	end do

        d0hi=1.0d40
        d0lo=1.0d40
        call spline(zrec,xrec,nz,d0lo,d0hi,dxrec)

        return
	end

C	===============================================================
	subroutine GET_INIT(z,x_H0,x_He0,x0)

C	Set up the initial conditions so it will work for general,
C	but not pathological choices of zstart
C	Initial ionization fraction using Saha for relevant species

	implicit none

	double precision OmegaT,HO,OmegaL,OmegaK
	double precision z,x0,rhs,x_H0,x_He0
	double precision Tnow,Nnow,z_eq
	double precision Lambda,H_frac
	double precision Lambda_He,Bfact,CK_He,CL_He
	double precision CB1,CDB,CR,CK,CL,CT,fHe,CB1_He1,CB1_He2,CDB_He,fu

	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He,fu
	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C	===============================================================

	if(z.gt.8000.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    x0 = 1.d0+2.d0*fHe

	else if(z.gt.3500.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He2/(Tnow*(1.d0+z)) ) / Nnow
	rhs = rhs*1.d0		!ratio of g's is 1 for He++ <-> He+
	x0 = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2
	1	+ 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )

	else if(z.gt.2000.d0)then

	x_H0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He1/(Tnow*(1.d0+z)) ) / Nnow
	rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs )
	1	- (rhs-1.d0))
	    x0 = x_He0
	    x_He0 = (x0 - 1.d0)/fHe

	else

	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1/(Tnow*(1.d0+z)) ) / Nnow
	    x_H0 = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )
	    x_He0 = 0.d0
	    x0 = x_H0

	end if

	return

	end

C	===============================================================
	subroutine ION(Ndim,z,Y,f)

	implicit none

	integer Ndim

	double precision z,x,n,n_He,Trad,Tmat,x_H,x_He
	double precision y(Ndim),f(Ndim)
	double precision C,k_B,h_P,m_e,m_H,sigma,a,Pi
	double precision Lambda,H_frac,Lambda_He
	double precision Tnow,HO,Nnow,z_eq,Hz,OmegaT,OmegaL,OmegaK
	double precision Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
	double precision timeTh,timeH
	double precision CB1,CDB,CR,CK,CL,CT,fHe,CB1_He1,CB1_He2,CDB_He,fu
	double precision Bfact,CK_He,CL_He
	double precision a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB

	common/Cfund/C,k_B,h_P,m_e,m_H,sigma,a,Pi
	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He,fu
	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C	===============================================================

c	the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen	
	a_PPB = 4.309d0
	b_PPB = -0.6166d0
	c_PPB = 0.6703d0
	d_PPB = 0.5300d0
c	the Verner and Ferland type fitting parameters for Helium
c	fixed to match those in the SSS papers, and now correct
	a_VF = 10.d0**(-16.744d0)
	b_VF = 0.711
	T_0 = 10.d0**(0.477121d0)	!3K
	T_1 = 10.d0**(5.114d0)

	x_H = y(1)
	x_He = y(2)
	x = x_H + fHe * x_He
	Tmat = y(3)

	n = Nnow * (1.d0+z)**3
	n_He = fHe * Nnow * (1.d0+z)**3
	Trad = Tnow * (1.d0+z)
	Hz = HO * dsqrt((1.d0+z)**4/(1.d0+z_eq)*OmegaT + OmegaT*(1.d0+z)**3
	1	+ OmegaK*(1.d0+z)**2 + OmegaL)

c	Get the radiative rates using PPQ fit (identical to Hummer's table)
	Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB
	1	/(1.d0+c_PPB*(Tmat/1.d4)**d_PPB)
	Rup = Rdown * (CR*Tmat)**(1.5d0)*dexp(-CDB/Tmat)

c	calculate He using a fit to a Verner & Ferland type formula
	sq_0 = dsqrt(Tmat/T_0)
	sq_1 = dsqrt(Tmat/T_1)
c	typo here corrected by Wayne Hu and Savita Gahlaut
	Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
	Rdown_He = rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
	Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*dexp(-CDB_He/Tmat)
	Rup_He = 4.d0*Rup_He	!statistical weights factor for HeI
c	Avoid overflow (pointed out by Jacques Roland)
	if((Bfact/Tmat).gt.680.d0)then
	  He_Boltz = dexp(680.d0)
	else
	  He_Boltz = dexp(Bfact/Tmat)
	end if
	K = CK/Hz		!Peebles coefficient K=lambda_a^3/8piH
	K_He = CK_He/Hz		!Peebles coefficient for Helium

c	Estimates of Thomson scattering time and Hubble time
	timeTh=(1.d0/(CT*Trad**4))*(1.d0+x+fHe)/x	!Thomson time
	timeH=2.d0/(3.d0*HO*(1.d0+z)**1.5)			!Hubble time

c	calculate the derivatives
c	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
c	(clunky, but seems to work)
	if (x_H.gt.0.99d0) then	!use Saha rate for Hydrogen
		f(1) = 0.d0
	else if (x_H.gt.0.98d0) then
		f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat))
	1	/(Hz*(1.d0+z))
	else
		f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat))
	1		*(1.d0 + K*Lambda*n*(1.d0-x_H)))
	2	/(Hz*(1.d0+z)*(1.d0/fu+K*Lambda*n*(1.d0-x)/fu
	3	+K*Rup*n*(1.d0-x)))
	end if
c	turn off the He once it is small
	if (x_He.lt.1.d-15) then
		f(2)=0.d0
	else
		f(2) = ((x*x_He*n*Rdown_He 
	1	- Rup_He*(1.d0-x_He)*dexp(-CL_He/Tmat))
	2		*(1. d0+ K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz))
	3	/(Hz*(1.d0+z)
	4	 * (1.d0 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))
	end if

c	follow the matter temperature once it has a chance of diverging
	if (timeTh.lt.H_frac*timeH) then
		f(3)=Tmat/(1.d0+z)	!Tmat follows Trad
	else
		f(3)= CT * (Trad**4) * x / (1.d0+x+fHe)
	1		* (Tmat-Trad) / (Hz*(1.d0+z)) + 2.d0*Tmat/(1.d0+z)
	end if

	return

	end
