#!/usr/bin/env python

from ctypes import *


class cmb(object):
	def __init__(self, shared_object="./libcmb.so.1.0"):
		self.libcmb = cdll.LoadLibrary(shared_object)
	def output(self, clts,cltt,cles,clet,clbt,clcs,clct,itflag,lmx):
		pass
	def thermo(self, tau, cs2b, opac, dxe):
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
c = cmb()
print c.indexx([3.,2.,1.])

