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
						
# Reynolds' stresses at building location
uu_b[:,0] = numpy.loadtxt(os.path.join(test0, 
           'postProcessing/sampleDict/'+timeStep[-1],'line_UPrime2Mean.xy'))[:,1]
vv_b[:,0] = numpy.loadtxt(os.path.join(test0, 
           'postProcessing/sampleDict/'+timeStep[-1],'line_UPrime2Mean.xy'))[:,4]
ww_b[:,0] = numpy.loadtxt(os.path.join(test0, 
           'postProcessing/sampleDict/'+timeStep[-1],'line_UPrime2Mean.xy'))[:,6]

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

#%% OPTIMIZATION LOOP
it = it + 1 

while (( (numpy.any(uu_b[:,it-1] < uu_exp_min) or 
          numpy.any(uu_b[:,it-1] > uu_exp_min) ) or

         (numpy.any(vv_b[:,it-1] < vv_exp_min) or 
          numpy.any(vv_b[:,it-1] > vv_exp_min) ) or
		  
		  
         (numpy.any(ww_b[:,it-1] < ww_exp_min) or 
          numpy.any(ww_b[:,it-1] > ww_exp_min) )) and
		  
          it < it_max):

	# test --------------------------------------------------------------------
	# case directory
	case = "/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/XCDF"+str(it)
	os.chdir(case)
	
	# Reynolds' stresses at building location
	[uu_b[:,it], vv_b[:,it], ww_b[:,it], Tu_b_mean[it], Tv_b_mean[it], Tw_b_mean[it]] = runOpenFOAM(case, Nproc, solver)
	
	# GRADIENT-BASED OPTIMIZATION ---------------------------------------------
	#uu_P = singleObjGrad(uu_b, uu_exp, uu_P, index_sample, relax_u, it)
	#vv_P = singleObjGrad(vv_b, vv_exp, vv_P, index_sample, relax_v, it)
	#ww_P = singleObjGrad(ww_b, ww_exp, ww_P, index_sample, relax_w, it)
	[vv_P, ww_P] = multiObjGrad_Rij(uu_b, uu_exp, vv_b, vv_exp, ww_b, ww_exp, vv_P, ww_P, index_probes, it)
	[Tv_0, Tw_0] = multiObjGrad_Ti(Tu_b_mean, 0.1482, Tv_b_mean, 0.0337, Tw_b_mean, 0.0469, Tv_0, Tw_0, it)
                      
	# PREPARE NEXT STEP -------------------------------------------------------
	# Create next step dir
	os.mkdir(os.path.join(case,'../XCDF' + str(it+1)))
	
	# Copy 0/constant/system directories
	shutil.copytree(os.path.join(test0,'0'), 
	                    os.path.join(case,'../XCDF' + str(it+1), '0'))
	shutil.copytree(os.path.join(test0,'constant'), 
	                    os.path.join(case,'../XCDF' + str(it+1), 'constant'))
	shutil.copytree(os.path.join(test0,'system'), 
	                    os.path.join(case,'../XCDF' + str(it+1), 'system'))
	
	# Bezier curves
	[uu_0, zuu_0] = bezier(uu_P[:,:,it+1], Nsample)
	[vv_0, zvv_0] = bezier(vv_P[:,:,it+1], Nsample)
	[ww_0, zww_0] = bezier(ww_P[:,:,it+1], Nsample)
	
	# Write input for new case    
	#uuBarInlet = numpy.array([zuu_0, uu_0]).T 
	vvBarInlet = numpy.array([zvv_0, vv_0]).T 
	wwBarInlet = numpy.array([zww_0, ww_0]).T
					 
	#numpy.savetxt(os.path.join(case,'../XCDF' + str(it+1), 
	#    'constant/uuBarInlet'), uuBarInlet, header='(', footer=')', comments='')
	numpy.savetxt(os.path.join(case,'../XCDF' + str(it+1), 
	    'constant/vvBarInlet'), vvBarInlet, header='(', footer=')', comments='')
	numpy.savetxt(os.path.join(case,'../XCDF' + str(it+1), 
	    'constant/wwBarInlet'), wwBarInlet, header='(', footer=')', comments='')
	
	inflowProp = ParsedParameterFile(os.path.join(case, "constant/inflowProperties"))
	#inflowProp["lagT_u"] = float(Tu_0[it])
	inflowProp["lagT_v"] = float(Tv_0[it+1])
	inflowProp["lagT_w"] = float(Tw_0[it+1])
	inflowProp.writeFile()
	
	# Update
	it = it + 1
