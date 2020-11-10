#------------------------------------------------------------------------------
#                      SIMPLE GRADIENT-BASED OPTIMIZATION
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# Load modules 
import numpy

def singleObjGrad(F, F_exp, x, index, relax, it):

	# Functions that given Reynolds' stresses or length-scales at 
	# building location, perform single objective gradient-based optimization
	# to compute cellZone parameters for the next step.
	#
	#     F        --> [...| F(z,it-1) | F(z,it) ] at building  
	#     F_exp    --> [...| F(z,it-1) | F(z,it) ] from experiment
	#     it       --> iteration number
	#     index    --> array containing indices corresponding to heights of
	#                    Bezier point
	#     relax    --> relaxation factor
	#     x        --> decision variable (Bezier points [Rij, z, it])
	
	Np = x.shape[0]
	
	# constraint
	for i in range(0,Np-2):
		
		# Newton method
		R  = F[index[i],it] - F_exp[index[i]]
		dR = F[index[i],it] - F[index[i],it-1]
			
		# Gradient of the unconstrained objective function
		gj = dR/(x[i+1,0,it] - x[i+1,0,it-1])

		# gradient of f
		df = R*gj
		
		# laplacian of f
		ddf = gj**2
		
		# Gauss-Newton method
		x[i+1,0,it+1] = x[i+1,0,it] - relax*df/ddf
		
		# check if R11, R22, R33 are positive
		if (x[i,0,it+1] < 0):
			x[i,0,it+1] = x[i,0,0]
			
	return x
