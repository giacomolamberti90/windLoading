#------------------------------------------------------------------------------
#                              OPENFOAM RUN
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# load modules
import os, numpy
from PyFoam.Applications.Decomposer           import Decomposer
from PyFoam.Applications.Runner               import Runner
from PyFoam.RunDictionary.SolutionDirectory   import SolutionDirectory
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from inputData                                import *
from timeScale        import lagrangianTimeScale

def runOpenFOAM(case, Nproc, solver):

	# Function that given number of processors, case directory and name of 
	# solver, decomposes, solves and checks if the simulation is converged 
	# and in that case returns the statistics at the building location
	#
	#     Rij_end --> Rij(z,it) at building  
	
	# convergence parameters
	tol     = 1e-3
	it      = 1
	it_max  = 2
	err     = 1
	
	# inizialization
	uu = numpy.zeros((Nprobes, it_max))
	vv = numpy.zeros((Nprobes, it_max))
	ww = numpy.zeros((Nprobes, it_max))
	
	uu_end = numpy.zeros((Nprobes, 1))
	vv_end = numpy.zeros((Nprobes, 1))
	ww_end = numpy.zeros((Nprobes, 1))
	
	Tu_end = numpy.zeros((Nprobes, 1))
	Tv_end = numpy.zeros((Nprobes, 1))
	Tw_end = numpy.zeros((Nprobes, 1))
	
	# time steps
	sol = SolutionDirectory(os.path.join(case,'postProcessing/probes/'))
	timeStep = []
	for t in sol:
		timeStep.append(t.baseName()); 

	# Reynolds stresses at building location
	# remove parentheses from file
	lines = open(os.path.join(case,'postProcessing/probes/',timeStep[-1],'UPrime2Mean')).readlines()[-1]
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
	uu_end = numpy.asarray(data_R[1::6])
	vv_end = numpy.asarray(data_R[4::6])
	ww_end = numpy.asarray(data_R[6::6])
	'''
	# integral time-scales
	data_L = []
	for t in timeStep[-5:]:
		
		# Lagrangian time-scales at building location
		# remove parentheses from file
		lines = open(os.path.join(case,'postProcessing/probes/',t,'U')).readlines()[Nprobes+2:]
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
		Tu_end[i] = lagrangianTimeScale(u, time-time[0])
		Tv_end[i] = lagrangianTimeScale(v, time-time[0])
		Tw_end[i] = lagrangianTimeScale(w, time-time[0])
	'''	
	# mean below 2m
	Tu_mean = 1;#numpy.mean(Tu_end[0:29])
	Tv_mean = 1;#numpy.mean(Tv_end[0:29])
	Tw_mean = 1;#numpy.mean(Tw_end[0:29])
	
	'''
	# Decompose 
	#Decomposer(args = [case , str(Nproc)])
	
	# CONVERGENCE LOOP --------------------------------------------------------
	while (err > tol and it < it_max):
		
		# Solve
		#Runner(args = ["--proc=%d" %Nproc, solver, "-case", case])
		
		# Final time
		te = SolutionDirectory(os.path.join(case, 'processor0')).getLast()
	
		# Reynolds' stresses at building location
		uu[:,it] = numpy.loadtxt(os.path.join(case, 
		       'postProcessing/sample/'+te, 'line_UPrime2Mean.xy'))[:,1]
		vv[:,it] = numpy.loadtxt(os.path.join(case, 
		       'postProcessing/sample/'+te, 'line_UPrime2Mean.xy'))[:,4]
		ww[:,it] = numpy.loadtxt(os.path.join(case, 
		       'postProcessing/sample/'+te, 'line_UPrime2Mean.xy'))[:,6]
		       		       
		# update   
		uu_end = uu[:,it]
		vv_end = vv[:,it]
		ww_end = ww[:,it]
		
		err_uu = numpy.amax(abs(uu[:,it] - uu[:,it-1])/abs(uu[:,it])) 
		err_vv = numpy.amax(abs(vv[:,it] - vv[:,it-1])/abs(vv[:,it])) 
		err_ww = numpy.amax(abs(ww[:,it] - ww[:,it-1])/abs(ww[:,it])) 
		err    = numpy.amax([err_uu, err_vv, err_ww, err_uv])  

		# change controlDict
		#controlDict = ParsedParameterFile("system/controlDict")
		#controlDict["endTime"] = controlDict["endTime"] + 2
		#controlDict.writeFile()
		
		# update
		it = it + 1	
	'''
	
	return (uu_end, vv_end, ww_end, Tu_mean, Tv_mean, Tw_mean)
	
	
