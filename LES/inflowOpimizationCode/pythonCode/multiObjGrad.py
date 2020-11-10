#------------------------------------------------------------------------------
#              MULTI-OBJECTIVE GRADIENT-BASED OPTIMIZATION
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# Load modules 
import numpy
from inputData  import relax_u, relax_v, relax_w, gamma_u, gamma_v, gamma_w

# Reynolds stresses optimization
def multiObjGrad_Rij(Fu, Fu_exp, Fv, Fv_exp, Fw, Fw_exp, x, y, index, it):

	# Functions that given Reynolds' stresses and/or lenght-scales at 
	# building location, perform multi-pbjective gradient-based optimization
	# to compute cellZone parameters for the next step.
	#
	#     Fi         --> [...| Fi(z,it-1) | Fi(z,it) ] at building  
	#     Fi_exp     --> [...| Fi(z,it-1) | Fi(z,it) ] from experiment
	#     it         --> iteration number
	#     index      --> array containing indices corresponding to heights of
	#                    Bezier point
	#     relax      --> relaxation factor
	#     x, y       --> decision variables (Bezier points [Rij, z, it])
	
	Np = x.shape[0]
	J  = numpy.zeros([2,1])
	H  = numpy.zeros([2,2])
	
	x[1,0,it+1] = x[1,0,1]
	y[1,0,it+1] = y[1,0,1]
   
	for i in range(1,Np-2):
		
		# Newton method
		Ru  = Fu[index[i],it] - Fu_exp[index[i]]
		dRu = Fu[index[i],it] - Fu[index[i],it-1]

		Rv  = Fv[index[i],it] - Fv_exp[index[i]]
		dRv = Fv[index[i],it] - Fv[index[i],it-1]

		Rw  = Fw[index[i],it] - Fw_exp[index[i]]
		dRw = Fw[index[i],it] - Fw[index[i],it-1]
			
		# Gradient of the single objective functions
		dFudx = dRu/(x[i+1,0,it] - x[i+1,0,it-1])
		dFudy = dRu/(y[i+1,0,it] - y[i+1,0,it-1])
		
		dFvdx = dRv/(x[i+1,0,it] - x[i+1,0,it-1])
		dFvdy = dRv/(y[i+1,0,it] - y[i+1,0,it-1])
		
		dFwdx = dRw/(x[i+1,0,it] - x[i+1,0,it-1])
		dFwdy = dRw/(y[i+1,0,it] - y[i+1,0,it-1])
		
		# gradient of f
		dfdx = gamma_u*Ru*dFudx + gamma_v*Rv*dFvdx + gamma_w*Rw*dFwdx
		dfdy = gamma_u*Ru*dFudy + gamma_v*Rv*dFvdy + gamma_w*Rw*dFwdy
   
		# second derivatives (Gauss-Newton)
		ddfddx = gamma_u*dFudx**2 + gamma_v*dFvdx**2 + gamma_w*dFwdx**2
		ddfddy = gamma_u*dFudy**2 + gamma_v*dFvdy**2 + gamma_w*dFwdy**2
		
		# Gauss-Newton method
		x[i+1,0,it+1] = x[i+1,0,it] - relax_v*dfdx/ddfddx
		y[i+1,0,it+1] = y[i+1,0,it] - relax_w*dfdy/ddfddy		
			
	return (x,y)

# integral time-scales optimization
def multiObjGrad_Ti(Fu, Fu_exp, Fv, Fv_exp, Fw, Fw_exp, x, y, it):

	# Functions that given Reynolds' stresses and/or lenght-scales at 
	# building location, perform multi-pbjective gradient-based optimization
	# to compute cellZone parameters for the next step.
	#
	#     Fi         --> [...| Fi(z,it-1) | Fi(z,it) ] at building  
	#     Fi_exp     --> [...| Fi(z,it-1) | Fi(z,it) ] from experiment
	#     it         --> iteration number
	#     index      --> array containing indices corresponding to heights of
	#                    Bezier point
	#     relax      --> relaxation factor
	#     x, y       --> decision variables (Bezier points [Rij, z, it])

	J  = numpy.zeros([2,1])
	H  = numpy.zeros([2,2])	
	
	# Newton method
	Ru  = Fu[it] - Fu_exp
	dRu = Fu[it] - Fu[it-1]

	Rv  = Fv[it] - Fv_exp
	dRv = Fv[it] - Fv[it-1]

	Rw  = Fw[it] - Fw_exp
	dRw = Fw[it] - Fw[it-1]
		
	# Gradient of the single objective functions
	dFudx = dRu/(x[it] - x[it-1])
	dFudy = dRu/(y[it] - y[it-1])
	
	dFvdx = dRv/(x[it] - x[it-1])
	dFvdy = dRv/(y[it] - y[it-1])
	
	dFwdx = dRw/(x[it] - x[it-1])
	dFwdy = dRw/(y[it] - y[it-1])
	
	# gradient of f
	dfdx = gamma_u*Ru*dFudx + gamma_v*Rv*dFvdx + gamma_w*Rw*dFwdx
	dfdy = gamma_u*Ru*dFudy + gamma_v*Rv*dFvdy + gamma_w*Rw*dFwdy
	
	# second derivatives
	ddfddx = gamma_u*dFudx**2 + gamma_v*dFvdx**2 + gamma_w*dFwdx**2
	ddfddy = gamma_u*dFudy**2 + gamma_v*dFvdy**2 + gamma_w*dFwdy**2
	
	# Gauss-Newton method
	x[it+1] = x[it] - relax_v*dfdx/ddfddx
	y[it+1] = y[it] - relax_w*dfdy/ddfddy
			
	return (x,y)
	
