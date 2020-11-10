#------------------------------------------------------------------------------
#                         LAGRANGIAN TIME-SCALES
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# load modules
import numpy
from numpy          import inf
from scipy.optimize import curve_fit

def func(x,a):
    return numpy.exp(-x/a)

def lagrangianTimeScale(x, time):
	
	N = len(x)
	
	# fluctuation components
	xf = x-numpy.mean(x)
	
	# autocovariance functions
	Rxx = numpy.correlate(xf, xf, "full")
	
	# autocorrelation functions
	RHOxx = Rxx[N-1:]/Rxx[N-1]
	
	# fit exponential function
	Tx = curve_fit(func, time, RHOxx, p0=1, bounds=(0,inf))[0]
		
	return (Tx)

