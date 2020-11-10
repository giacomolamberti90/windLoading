#------------------------------------------------------------------------------
#                              POST-PROCESSING
#------------------------------------------------------------------------------

# load modules
import os, numpy, matplotlib
import matplotlib.pyplot as plt

from   PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory

Nsample = 49
it      = 2

# height
z0 = 0.0032
z  = numpy.linspace(0.04, 1.96, Nsample)

# initialization
U  = numpy.zeros((Nsample, 5, it))
uu = numpy.zeros((Nsample, 5, it))
vv = numpy.zeros((Nsample, 5, it))
ww = numpy.zeros((Nsample, 5, it))
uv = numpy.zeros((Nsample, 5, it))

uu_0     = numpy.zeros((Nsample, it))
vv_0     = numpy.zeros((Nsample, it))
ww_0     = numpy.zeros((Nsample, it))
uv_0     = numpy.zeros((Nsample, it))

uu_end   = numpy.zeros((Nsample, it))
uu_end_1 = numpy.zeros((Nsample, it))

leg = []

#%% GOAL: EXPERIMENTAL DATA AT BUILDING LOCATION 
exp = "/home/giacomol/Desktop/optimization/exp"

# mean
U_exp = numpy.genfromtxt(os.path.join(exp,'UInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp = numpy.genfromtxt(os.path.join(exp,'uuBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp = numpy.genfromtxt(os.path.join(exp,'vvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp = numpy.genfromtxt(os.path.join(exp,'wwBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp = numpy.genfromtxt(os.path.join(exp,'uvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

# mean - 1.96*std
U_exp_min = numpy.genfromtxt(os.path.join(exp,'UInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp_min = numpy.genfromtxt(os.path.join(exp,'uuBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp_min = numpy.genfromtxt(os.path.join(exp,'vvBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp_min = numpy.genfromtxt(os.path.join(exp,'wwBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp_min = numpy.genfromtxt(os.path.join(exp,'uvBarInlet_min'), 
                                usecols=1, skip_header=1, skip_footer=1)

# mean + 1.96*std
U_exp_max = numpy.genfromtxt(os.path.join(exp,'UInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
uu_exp_max = numpy.genfromtxt(os.path.join(exp,'uuBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
vv_exp_max = numpy.genfromtxt(os.path.join(exp,'vvBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
ww_exp_max = numpy.genfromtxt(os.path.join(exp,'wwBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)
uv_exp_max = numpy.genfromtxt(os.path.join(exp,'uvBarInlet_max'), 
                                usecols=1, skip_header=1, skip_footer=1)


te = ['125','130']

#%% Statistics at building location
for j in  range(0,it):
	
	case = '/home/giacomol/Desktop/optimization/test0'
	#te   = SolutionDirectory(os.path.join(case)).getLast()
	legj = 'step' + str(j)
	leg.append(legj)

	# cellZone
	uu_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uuBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	vv_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/vvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	ww_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/wwBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	uv_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	# Building location 
	for i in range(0,5):		
		
	    U[:,i,j] = numpy.loadtxt(os.path.join(case, 
	    'postProcessing/sampleDict/'+te[j], 'line'+str(i)+'_UMean.xy'))[:,1]
		
	    uu[:,i,j] = numpy.loadtxt(os.path.join(case, 
	    'postProcessing/sampleDict/'+te[j], 'line'+str(i)+'_UPrime2Mean.xy'))[:,1]
		
	    vv[:,i,j] = numpy.loadtxt(os.path.join(case, 
	    'postProcessing/sampleDict/'+te[j], 'line'+str(i)+'_UPrime2Mean.xy'))[:,4]	
		
	    ww[:,i,j] = numpy.loadtxt(os.path.join(case, 
	    'postProcessing/sampleDict/'+te[j], 'line'+str(i)+'_UPrime2Mean.xy'))[:,6]	
	
	    uv[:,i,j] = numpy.loadtxt(os.path.join(case, 
	    'postProcessing/sampleDict/'+te[j], 'line'+str(i)+'_UPrime2Mean.xy'))[:,2]				
	
# spanwise average
U_mean  = numpy.average(U,1)
uu_mean = numpy.average(uu,1)
vv_mean = numpy.average(vv,1)
ww_mean = numpy.average(ww,1)
uv_mean = numpy.average(uv,1)

#%% pressure fluctuations

pp = numpy.zeros((Nsample, 5))

# streamwise
for i in range(0,5):		
	
    pp[:,i] = numpy.loadtxt(os.path.join(case, 
    'postProcessing/sampleDict/'+te[1], 'line'+str(i)+'_pPrime2Mean.xy'))[:,1]
	
# spanwise average
pp_mean  = numpy.average(pp,1)

#%% FIGURES

fig = '/home/giacomol/Desktop/Research/slides/optimization/figures/'
matplotlib.rcParams.update({'font.size': 14})

leg2 = leg
leg2.append('exp')
leg3 = leg2
leg3.append('Inlet')

#PRESSURE FLUCTUATIONS ALONG WIND TUNNEL ----------------------------------------
plt.plot(numpy.linspace(-5,5,Nsample), pp_mean)
plt.xlabel('x')
plt.ylabel('pp', rotation=0)
#plt.savefig(os.path.join(fig, 'U.png'))
plt.show()

# BUILDING --------------------------------------------------------------------
plt.plot(U_mean, z)
plt.plot(U_exp, z+z0, 'ok', alpha=.2)	
plt.fill_betweenx(z+z0, U_exp_min, U_exp_max, color='k', alpha=.2)
plt.xlabel('U')
plt.ylabel('z', rotation=0)
plt.legend(leg2)
#plt.savefig(os.path.join(fig, 'U.png'))
plt.show()

plt.plot(uu_mean, z)
plt.plot(uu_exp, z+z0, 'ok', alpha=.2)	
plt.fill_betweenx(z+z0, uu_exp_min, uu_exp_max, color='k', alpha=.2)
#plt.plot(uu_0, z)	
plt.xlabel('uu')
plt.ylabel('z', rotation=0)
plt.legend(leg2)
#plt.savefig(os.path.join(fig, 'uu.png'))
plt.show()

plt.plot(vv_mean, z)	
plt.plot(vv_exp, z+z0, 'ok', alpha=.2)	
plt.fill_betweenx(z+z0, vv_exp_min, vv_exp_max, color='k', alpha=.2)
#plt.plot(vv_0, z)	
plt.xlabel('vv')
plt.ylabel('z', rotation=0)
#plt.savefig(os.path.join(fig, 'vv.png'))
plt.show()

plt.plot(ww_mean, z)	
plt.plot(ww_exp, z+z0, 'ok', alpha=.2)	
plt.fill_betweenx(z+z0, ww_exp_min, ww_exp_max, color='k', alpha=.2)
#plt.plot(ww_0, z)	
plt.xlabel('ww')
plt.ylabel('z', rotation=0)
#plt.savefig(os.path.join(fig, 'ww.png'))
plt.show()

plt.plot(uv_mean, z)	
plt.plot(uv_exp, z+z0, 'ok', alpha=.2)	
plt.fill_betweenx(z+z0, uv_exp_min, uv_exp_max, color='k', alpha=.2)
#plt.plot(uv_0, z)	
plt.xlabel('uv')
plt.ylabel('z', rotation=0)
#plt.savefig(os.path.join(fig, 'uv.png'))
plt.show()

'''
# CELLZONE --------------------------------------------------------------------
plt.plot(uu_0, z)	
plt.xlabel('uuBarInlet')
plt.ylabel('z', rotation=0)
plt.legend(leg)
#plt.savefig(os.path.join(fig, 'uuBarInlet.eps'))
plt.show()

plt.plot(vv_0, z)	
plt.xlabel('vvBarInlet')
plt.ylabel('z', rotation=0)
#plt.savefig(os.path.join(fig, 'vvBarInlet.eps'))
plt.show()

plt.plot(ww_0, z)	
plt.xlabel('wwBarInlet')
plt.ylabel('z', rotation=0)
#plt.savefig(os.path.join(fig, 'wwBarInlet.eps'))
plt.show()
'''

