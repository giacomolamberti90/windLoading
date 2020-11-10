#------------------------------------------------------------------------------
#                           USER-DEFINED INPUTS
#------------------------------------------------------------------------------

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:13:41 2017

@author: Giacomo Lamberti
"""

# Load modules 
import os, numpy

#%% GOAL: EXPERIMENTAL DATA AT BUILDING LOCATION ------------------------------

# directory with mean and 95 % CI of experimental mean velocity 
# and Reynolds' stresses
exp = "/home/storage/LES/optimization/highRiseDomain/exp/"

zsample = numpy.genfromtxt(os.path.join(exp,'UInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)
zprobes = numpy.genfromtxt(os.path.join(exp,'TuInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)

# mean
U_exp  = numpy.genfromtxt(os.path.join(exp,'UInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp = numpy.genfromtxt(os.path.join(exp,'uuBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp = numpy.genfromtxt(os.path.join(exp,'vvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp = numpy.genfromtxt(os.path.join(exp,'wwBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp = numpy.genfromtxt(os.path.join(exp,'uvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

Tu_exp = numpy.genfromtxt(os.path.join(exp,'TuInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tv_exp = numpy.genfromtxt(os.path.join(exp,'TvInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tw_exp = numpy.genfromtxt(os.path.join(exp,'TwInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

xLu_exp = numpy.genfromtxt(os.path.join(exp,'xLuInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLv_exp = numpy.genfromtxt(os.path.join(exp,'xLvInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLw_exp = numpy.genfromtxt(os.path.join(exp,'xLwInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

# mean - 1.96*std
U_exp_min  = numpy.genfromtxt(os.path.join(exp,'UInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp_min = numpy.genfromtxt(os.path.join(exp,'uuBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp_min = numpy.genfromtxt(os.path.join(exp,'vvBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp_min = numpy.genfromtxt(os.path.join(exp,'wwBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp_min = numpy.genfromtxt(os.path.join(exp,'uvBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)

Tu_exp_min = numpy.genfromtxt(os.path.join(exp,'TuInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tv_exp_min = numpy.genfromtxt(os.path.join(exp,'TvInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tw_exp_min = numpy.genfromtxt(os.path.join(exp,'TwInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)

xLu_exp_min = numpy.genfromtxt(os.path.join(exp,'xLuInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLv_exp_min = numpy.genfromtxt(os.path.join(exp,'xLvInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLw_exp_min = numpy.genfromtxt(os.path.join(exp,'xLwInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)

# mean + 1.96*std
U_exp_max  = numpy.genfromtxt(os.path.join(exp,'UInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp_max = numpy.genfromtxt(os.path.join(exp,'uuBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp_max = numpy.genfromtxt(os.path.join(exp,'vvBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp_max = numpy.genfromtxt(os.path.join(exp,'wwBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp_max = numpy.genfromtxt(os.path.join(exp,'uvBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)

Tu_exp_max = numpy.genfromtxt(os.path.join(exp,'TuInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tv_exp_max = numpy.genfromtxt(os.path.join(exp,'TvInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
Tw_exp_max = numpy.genfromtxt(os.path.join(exp,'TwInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)

xLu_exp_max = numpy.genfromtxt(os.path.join(exp,'xLuInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLv_exp_max = numpy.genfromtxt(os.path.join(exp,'xLvInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
xLw_exp_max = numpy.genfromtxt(os.path.join(exp,'xLwInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)

# power spectrum at 1m
Euu_exp = numpy.genfromtxt(os.path.join(exp,'Euu_exp'), usecols=1)

# turbulence kinetic energy
k_exp     = 0.5 * (uu_exp + vv_exp + ww_exp)
k_exp_min = 0.5 * (uu_exp_min + vv_exp_min + ww_exp_min)
k_exp_max = 0.5 * (uu_exp_max + vv_exp_max + ww_exp_max)

#%% OpenFOAM PARAMETERS -------------------------------------------------------

# number of processors
Nproc = 8

# LES solver
solver = "myXCDF"

# number of points 
Nsample = 56
Nprobes = 50

#%% OPTIMIZATION PARAMETERS ---------------------------------------------------

# multi-objective weights
gamma_u = 0.5
gamma_v = 0.2
gamma_w = 0.3

# tolerance
tol = 1e-3

# max number of iterations
it_max = 4

# relaxation factors
relax_u = 1.0
relax_v = 1.0
relax_w = 1.0

# directories of initial simulations (2 because of gradient)
test0 = "/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/XCDF0" 
test1 = "/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/XCDF1" 

# directory of first optimization step
test2 = "/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/XCDF2" 

#%% BEZIER PARAMETRIZATION PARAMETERS -----------------------------------------
# high-rise domain
# Number of points
Np = 8

# x-coordinate is fixed while y is optimized!
P_x = [0.04, 0.08, 0.30, 0.50, 1.00, 2.00, 3.00, 3.96]

# indices of P_x 
index_probes = [1, 6, 11, 19, 29, 41]

# Bezier points for uu (initial simulations)
P0_uu = [0.8,   1.5,  1.0,  0.8,  0.3,  0.1, 0.05,  0.1]
P1_uu = [0.8,   3.0,  2.0,  1.6,  0.6,  0.2,  0.1,  0.1]

# Bezier points for vv (initial simulations)
P0_vv = [0.1,   1.0,  0.6,  0.4,  0.3,  0.2, 0.05,  0.1]
P1_vv = [0.1,   2.0,  1.2,  0.8,  0.6,  0.4,  0.1,  0.1]

# Bezier points for ww (initial simulations)
P0_ww = [0.3,   1.2,  0.8,  0.6,  0.4,  0.2,  0.1,  0.1]
P1_ww = [0.3,   2.4,  1.6,  1.2,  0.8,  0.45,  0.2,  0.1]
