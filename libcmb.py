#!/usr/bin/env python

from ctypes import *
from struct import *


class cmb(object):
	def __init__(self, libcmb="./libcmb.so.1.0", libjlgen="./libjlgen.so.1.0"):
		self.libcmb = cdll.LoadLibrary(libcmb)
		self.libjlgen = cdll.LoadLibrary(libjlgen)
		
		# Set default values
		self.ict = 0
		self.lmoin = 1500
		self.akmax0 = 3000
		self.l0max = 5300 # Use the same value as in cmbfast.inc

		
		
	def __setattr__(self, name, value):
		"""Customize attribute asignment to guarantee sane values. 
		"""
		if name == "ict":
			if value in [0,1,2]:
				object.__setattr__(self, name, value)
			else:
				raise ValueError("Should be 0 (CMB), 1 (transfer functions) or 2 (both)")
		elif name == "lmoin":
			object.__setattr__(self, name, value)
		elif name == "akmax0":
			object.__setattr__(self, name, value)
		else: # Default
			object.__setattr__(self, name, value)
	
	def _gen_fortran_block_(self, fmt, v):
		""" A helper function to generate a fortran block when writing to file.
		"""
		header = pack("P", calcsize(fmt))
		return header + pack(fmt, v) + header
	
	def jlgen(self, lmoin ,kmax0, filename=False):
		if lmoin + 300 > self.l0max:
			raise ValueError("lmoin should be less than %s." % (self.l0max-299))
		self.libjlgen.initlval_(byref(c_int(lmoin)))
		if filename:
			f = open(filename, "wb")
		else: # Use auto filename
			f = open("jl%dx%d.dat" % (lmoin, kmax0), "wb")
		
		# Write lmo ( = lmoin + 300 )
		f.write(self._gen_fortran_block_("L", lmoin + 300))
		# write kmax0
		f.write(self._gen_fortran_block_("L", kmax0))

		f.close()
		
	
	# The subroutines.F functions and subroutines. 
	def output(self, clts,cltt,cles,clet,clbt,clcs,clct,itflag,lmx):
		pass
	def COBEnormalize(self, clts,cltt,cles,clet,clbt,clcs,clct,clkk,cltk):
		pass
	def initlval(self, lmoin,akmax0in):
		pass
	def thermo(self, tau, cs2b, opac, dxe):
		pass
	def ionize(self, tempb,a,adot,dtau,xe):
		pass
	def ionhe(self, tempb, a, x0, x1, x2):
		pass
	def nu1(self, a,rhonu,pnu):
		pass
	def initnu1(self, amnu):
		pass
	def ninu1(self, a,rhonu,pnu,amnu):
		pass
	def nu2(self, a,drhonu,fnu,dpnu,shearnu,psi0,psi1,psi2):
		pass
	def nuder(self, a,adotoa,rhonu,rhonudot,shearnudot,psi2,psi2dot):
		pass
	def dsoundda(self, a):
		pass
	def reiopar(self, optdlss,zri,zristp,rif):
		pass
	def recint(self, a,xe):
		pass
	def splder(self, y,dy,n):
		pass
	def rombint(self, f,a,b,tol):
		pass
	def spline(self, x,y,yp1,ypn):
		if len(x) != len(y):
			raise ValueError("x and y should be of the same length.")
		l = len(x)
		fn = c_int(l)
		fx = (c_double * l)()
		fy = (c_double * l)()
		for i in xrange(l):
			fx[i] = x[i]
			fy[i] = y[i]
		fn = c_int(l)
		fyp1 = c_double(yp1)
		fypn = c_double(ypn)
		fy2 = (c_double * l)()
		self.libcmb.spline_(fx, fy, byref(fn), byref(fyp1), byref(fypn), fy2)
		return list(fy2)

	def splint(self, y,z,n):
		pass
	def output_power(self, ntf,amnu):
		pass
	def outtransf(self, n,y,curv,ik,itf):
		pass
	def setuptransf(self, ntf):
		pass
	def out(self, a,b,c,n,file):
		pass
	def indexx(self, arr):
		l = len(arr)
		fn = c_int(l)
		farr = (c_double * l)()
		for i in xrange(l):
			farr[i] = arr[i]
		findx = (c_int * l)()
		self.libcmb.indexx_(byref(fn), farr, findx)
		return list(findx)
	def dynrho(self, a):
		pass
	def wdyn_func(self, a):
		pass
	def dwda_func(self, a):
		pass
	def dyn_phi(self,a,hdot,grho,gpres,phi,psi,dphi,dpsi):
		pass
	def dyn_nrg(self,a,grho,gpres,delphi,delpsi,dgrho,dgprs,dgtheta):
		pass
	def readtable(self):
		pass
	

if __name__ == "__main__":
	c = cmb()
	print c.jlgen(1500, 3000)

