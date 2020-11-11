#------------------------------------------------------------------------------
#                              OPTIMIZATION 
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# Load modules 
import os, shutil, numpy, numpy.matlib

from singleObjGrad    import singleObjGrad
from multiObjGrad     import *
from runOpenFOAM      import runOpenFOAM
from bezier           import bezier
from inputData        import *
from timeScale        import lagrangianTimeScale

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.RunDictionary.SolutionDirectory   import SolutionDirectory

#%% PREPARE INITIAL CASES

# Inizialization
# Mean velocity at building location
U_b = numpy.zeros((Nprobes, it_max))

# Reynolds' stresses at building location
uu_b = numpy.zeros((Nprobes, it_max))
vv_b = numpy.zeros((Nprobes, it_max))
ww_b = numpy.zeros((Nprobes, it_max))

# Lagrangian time-scales at building location
Tu_b = numpy.zeros((Nprobes, it_max))
Tv_b = numpy.zeros((Nprobes, it_max))
Tw_b = numpy.zeros((Nprobes, it_max))

# Lagrangian time-scales at building location mean below 2 meters
Tu_b_mean = numpy.zeros((it_max, 1))
Tv_b_mean = numpy.zeros((it_max, 1))
Tw_b_mean = numpy.zeros((it_max, 1))

# Bezier points to parametrize Reynolds' stresses at cellZone
uu_P = numpy.zeros((Np, 2, it_max))
vv_P = numpy.zeros((Np, 2, it_max))
ww_P = numpy.zeros((Np, 2, it_max))

# Lagrangian time-scales at cellZone
Tu_0 = numpy.zeros((it_max, 1))
Tv_0 = numpy.zeros((it_max, 1))
Tw_0 = numpy.zeros((it_max, 1))

#%% test0 (input from experiment) ----------------------------------------------
# BEZIER POINTS
#uu
uu_P[:,0,0]  = P0_uu      # optimize
uu_P[0,0,:]  = P0_uu[0]   # first point fixed
uu_P[-1,0,:] = P0_uu[-1]  # last point fixed
uu_P[:,1,:]  = numpy.matlib.repmat(P_x, it_max, 1).T
#vv
vv_P[:,0,0]  = P0_vv 
vv_P[0,0,:]  = P0_vv[0]
vv_P[-1,0,:] = P0_vv[-1]
vv_P[:,1,:]  = uu_P[:,1,:]
#ww
ww_P[:,0,0]  = P0_ww
ww_P[0,0,:]  = P0_ww[0]
ww_P[-1,0,:] = P0_ww[-1]
ww_P[:,1,:]  = uu_P[:,1,:]

# time steps
sol = SolutionDirectory(os.path.join(test0, 'postProcessing/probes/'))
timeStep = []
for t in sol:
	timeStep.append(t.baseName()); 

# Reynolds stresses at building location
# remove parentheses from file
lines = open(os.path.join(test0,'postProcessing/probes/',timeStep[-1],'UPrime2Mean')).readlines()[-1]
temp  = str(lines)
temp  = temp.replace('(','')	
temp  = temp.replace(')','')	

data_R = []
# read floats
for elem in temp.split():
        try:
            data_R.append(float(elem))
        except ValueError:
            pass

# Reynolds stresses
uu_b[:,0] = numpy.asarray(data_R[1::6])
vv_b[:,0] = numpy.asarray(data_R[4::6])
ww_b[:,0] = numpy.asarray(data_R[6::6])

# integral length-scales
Tv_0[0] = 0.03
Tw_0[0] = 0.04

data_L = []
for t in timeStep[-1:]:
	
	# Lagrangian time-scales at building location
	# remove parentheses from file
	lines = open(os.path.join(test0,'postProcessing/probes/',t,'U')).readlines()[Nprobes+2:]
	temp  = str(lines)
	temp  = temp.replace('(','')	
	temp  = temp.replace(')','')
	temp  = temp.replace('\\n',' ')	
	
	# read floats
	for elem in temp.split():
	        try:
	            data_L.append(float(elem))
	        except ValueError:
	            pass

# time
time = numpy.asarray(data_L[0::1+3*Nprobes])

for i in range(0,Nprobes):

	# extract data
	u = numpy.asarray(data_L[1+3*i::1+3*Nprobes])
	v = numpy.asarray(data_L[2+3*i::1+3*Nprobes])
	w = numpy.asarray(data_L[3+3*i::1+3*Nprobes])
	
	# time-scales at building location
	Tu_b[i,0] = lagrangianTimeScale(u, time-time[0])
	Tv_b[i,0] = lagrangianTimeScale(v, time-time[0])
	Tw_b[i,0] = lagrangianTimeScale(w, time-time[0])

# mean below 2m	
Tu_b_mean[0] = 0.1288
Tv_b_mean[0] = 0.0539
Tw_b_mean[0] = 0.0643
		 
#%% test1 (double input) -------------------------------------------------------
# BEZIER POINTS
uu_P[:,0,1] = P1_uu
vv_P[:,0,1] = P1_vv 
ww_P[:,0,1] = P1_ww

# time steps
sol = SolutionDirectory(os.path.join(test1, 'postProcessing/probes/'))
timeStep = []
for t in sol:
	timeStep.append(t.baseName()); 

# Reynolds stresses at building location
# remove parentheses from file
lines = open(os.path.join(test1,'postProcessing/probes/',timeStep[-1],'UPrime2Mean')).readlines()[-1]
temp  = str(lines)
temp  = temp.replace('(','')	
temp  = temp.replace(')','')	

data_R = []
# read floats
for elem in temp.split():
        try:
            data_R.append(float(elem))
        except ValueError:
			
            pass
# Reynolds stresses
uu_b[:,1] = numpy.asarray(data_R[1::6])
vv_b[:,1] = numpy.asarray(data_R[4::6])
ww_b[:,1] = numpy.asarray(data_R[6::6])
	
# integral length-scales
Tv_0[1] = 0.06
Tw_0[1] = 0.08
	
data_L = []
for t in timeStep[-1:]:
	
	# Lagrangian time-scales at building location
	# remove parentheses from file
	lines = open(os.path.join(test1,'postProcessing/probes/',t,'U')).readlines()[Nprobes+2:]
	temp  = str(lines)
	temp  = temp.replace('(','')	
	temp  = temp.replace(')','')	
	temp  = temp.replace('\\n',' ')	
	
	# read floats
	for elem in temp.split():
	        try:
	            data_L.append(float(elem))
	        except ValueError:
	            pass

# time
time = numpy.asarray(data_L[0::1+3*Nprobes])

for i in range(0,Nprobes):

	# extract data
	u = numpy.asarray(data_L[1+3*i::1+3*Nprobes])
	v = numpy.asarray(data_L[2+3*i::1+3*Nprobes])
	w = numpy.asarray(data_L[3+3*i::1+3*Nprobes])
	
	# time-scales at building location
	Tu_b[i,1] = lagrangianTimeScale(u, time-time[0])
	Tv_b[i,1] = lagrangianTimeScale(v, time-time[0])
	Tw_b[i,1] = lagrangianTimeScale(w, time-time[0])

# mean below 2m
Tu_b_mean[1] = 0.1570
Tv_b_mean[1] = 0.0668
Tw_b_mean[1] = 0.0935
		 
#%% test2 (first optimization step) -------------------------------------------
it = 1
		  
# compute input points
[vv_P, ww_P] = multiObjGrad_Rij(uu_b, uu_exp, vv_b, vv_exp, ww_b, ww_exp, vv_P, ww_P, index_probes, it)
[Tv_0, Tw_0] = multiObjGrad_Ti(Tu_b_mean, 0.1482, Tv_b_mean, 0.0337, Tw_b_mean, 0.0469, Tv_0, Tw_0, it)

# Bezier curves
[uu_0, zuu_0] = bezier(uu_P[:,:,2], Nsample)
[vv_0, zvv_0] = bezier(vv_P[:,:,2], Nsample)
[ww_0, zww_0] = bezier(ww_P[:,:,2], Nsample)

# Write input for test2		
numpy.savetxt(os.path.join(test2, 'constant/vvBarInlet'), 
	     numpy.array([zvv_0, vv_0]).T, header='(', footer=')', comments='')
numpy.savetxt(os.path.join(test2, 'constant/wwBarInlet'), 
	     numpy.array([zww_0, ww_0]).T, header='(', footer=')', comments='')


inflowProp = ParsedParameterFile(os.path.join(test2, "constant/inflowProperties"))
inflowProp["lagT_v"] = float(Tv_0[it+1])
inflowProp["lagT_w"] = float(Tw_0[it+1])
inflowProp.writeFile()
