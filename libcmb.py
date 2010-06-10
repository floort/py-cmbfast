#!/usr/bin/env python

from ctypes import *


class cmb(object):
	def __init__(self, shared_object="./libcmb.so.1.0"):
		self.libcmb = cdll.LoadLibrary(shared_object)
		
		# Set default values
		self.ict = 0
		self.lmoin = 1500
		self.akmax0 = 3000
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
	def spline(self, x,y,n,yp1,ypn,y2):
		pass
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
	print c.indexx([3.,2.,1.])

