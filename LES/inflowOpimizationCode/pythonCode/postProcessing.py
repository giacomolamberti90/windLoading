#------------------------------------------------------------------------------
#                              POST-PROCESSING
#------------------------------------------------------------------------------

# load modules
import os, numpy, matplotlib
import matplotlib.pyplot as plt
import scipy.io as spio

from timeScale                              import lagrangianTimeScale
from inputData                              import *
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
from scipy.interpolate                      import interp1d

it = 4
Tf = 20
dt = 0.0015
N  = int(Tf/dt)

# initialization
U  = numpy.zeros((Nprobes, it-1))
uu = numpy.zeros((Nprobes, it-1))
vv = numpy.zeros((Nprobes, it-1))
ww = numpy.zeros((Nprobes, it-1))
uv = numpy.zeros((Nprobes, it-1))

z_U  = numpy.zeros((Nsample, it))
z_uu = numpy.zeros((Nsample, it))
z_vv = numpy.zeros((Nsample, it))
z_ww = numpy.zeros((Nsample, it))
z_uv = numpy.zeros((Nsample, it))

U_0  = numpy.zeros((Nsample, it))
uu_0 = numpy.zeros((Nsample, it))
vv_0 = numpy.zeros((Nsample, it))
ww_0 = numpy.zeros((Nsample, it))
uv_0 = numpy.zeros((Nsample, it))

Tu = numpy.zeros((Nprobes, it-1))
Tv = numpy.zeros((Nprobes, it-1))
Tw = numpy.zeros((Nprobes, it-1))

xLu = numpy.zeros((Nprobes, it-1))
xLv = numpy.zeros((Nprobes, it-1))
xLw = numpy.zeros((Nprobes, it-1))

freq = numpy.zeros((N, it-1))
Eu   = numpy.zeros((N, it-1))
Ev   = numpy.zeros((N, it-1))
Ew   = numpy.zeros((N, it-1))

leg = []

#%% Statistics at cellZone location
for j in  range(0,it-1):
	
	case = os.path.join('/home/storage/LES/optimization/highRiseDomain/XCDF' +str(j)) 
	#case = os.path.join('/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/TKE.' +str(j))
	
	te   = SolutionDirectory(os.path.join(case)).getLast()
	legj = 'step' + str(j)
	leg.append(legj)
	
	# cellZone
	z_U[:,j] = numpy.genfromtxt(os.path.join(case,'constant/UInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)
	U_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/UInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	# cellZone
	z_uu[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uuBarInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)
	uu_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uuBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)
	
	z_vv[:,j] = numpy.genfromtxt(os.path.join(case,'constant/vvBarInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)
	vv_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/vvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

	z_ww[:,j] = numpy.genfromtxt(os.path.join(case,'constant/wwBarInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)	
	ww_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/wwBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

	z_uv[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uvBarInlet'), 
                                usecols=0, skip_header=1, skip_footer=1)	
	uv_0[:,j] = numpy.genfromtxt(os.path.join(case,'constant/uvBarInlet'), 
                                usecols=1, skip_header=1, skip_footer=1)

#%% Statistics at building location
for j in  range(0,it-1):
	
	case = os.path.join('/home/storage/LES/optimization/highRiseDomain/XCDF' +str(j)) 
	#case = os.path.join('/home/giacomol/Desktop/Research/windLoading/LES/optimization/highRiseDomain/TKE.' +str(j))
	
	# time steps
	sol = SolutionDirectory(os.path.join(case,'postProcessing/probes/'))
	timeStep = []
	for t in sol:
		timeStep.append(t.baseName()); 
				 
	# Reynolds stresses at building location
	# remove parentheses from file
	lines = open(os.path.join(case,'postProcessing/probes/',timeStep[-1],'UMean')).readlines()[-1]
	temp  = str(lines)
	temp  = temp.replace('(','')	
	temp  = temp.replace(')','')	
	
	data_U = []
	# read floats
	for elem in temp.split():
	        try:
	            data_U.append(float(elem))
	        except ValueError:
	            pass
	
	# Reynolds stresses
	U[:,j] = numpy.asarray(data_U[1::3])
	
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
	uu[:,j] = numpy.asarray(data_R[1::6])
	vv[:,j] = numpy.asarray(data_R[4::6])
	ww[:,j] = numpy.asarray(data_R[6::6])
	uv[:,j] = numpy.asarray(data_R[2::6])
	
	data_L = []
	for t in timeStep[-numpy.minimum(len(timeStep)-2,5):]:
		
		# Lagrangian time-scales at building location
		# remove parentheses from file
		lines = open(os.path.join(case,'postProcessing/probes/', t,'U')).readlines()[Nprobes+2:]
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
	print(time[-1]-time[0])
	
	for i in range(0,Nprobes):
	
		# extract data
		u = numpy.asarray(data_L[1+3*i::1+3*Nprobes])
		v = numpy.asarray(data_L[2+3*i::1+3*Nprobes])
		w = numpy.asarray(data_L[3+3*i::1+3*Nprobes])
		
		# time-scales at building location
		Tu[i,j] = lagrangianTimeScale(u, time-time[0])
		Tv[i,j] = lagrangianTimeScale(v, time-time[0])
		Tw[i,j] = lagrangianTimeScale(w, time-time[0])
	
	# integral length-scales:
	xLu[:,j] = Tu[:,j]*U[:,j]
	xLv[:,j] = Tv[:,j]*U[:,j]
	xLw[:,j] = Tw[:,j]*U[:,j]
	
# turbulence intensities
Iu_exp     = numpy.sqrt(uu_exp)/U_exp
Iu_exp_min = numpy.sqrt(uu_exp_min)/U_exp_max
Iu_exp_max = numpy.sqrt(uu_exp_max)/U_exp_min
					   
Iv_exp     = numpy.sqrt(vv_exp)/U_exp
Iv_exp_min = numpy.sqrt(vv_exp_min)/U_exp_max
Iv_exp_max = numpy.sqrt(vv_exp_max)/U_exp_min
					   
Iw_exp     = numpy.sqrt(ww_exp)/U_exp
Iw_exp_min = numpy.sqrt(ww_exp_min)/U_exp_max
Iw_exp_max = numpy.sqrt(ww_exp_max)/U_exp_min
					   
Iuv_exp     = uv_exp/U_exp**2
Iuv_exp_min = uv_exp_min/U_exp_min**2
Iuv_exp_max = uv_exp_max/U_exp_max**2

k         = 0.5*(uu + vv + ww) 
k_exp     = 0.5*(uu_exp + vv_exp + ww_exp)
k_exp_min = 0.5*(uu_exp_min + vv_exp_min + ww_exp_min)
k_exp_max = 0.5*(uu_exp_max + vv_exp_max + ww_exp_max)
					   
Iu  = numpy.sqrt(uu)/U
Iv  = numpy.sqrt(vv)/U
Iw  = numpy.sqrt(ww)/U
Iuv = uv/U**2
			   
#%% FIGURES

ustar = 0.49
H     = 2

mat = spio.loadmat('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/data/ycoord.mat', squeeze_me=True)
zexp = mat['ZHW']

mat = spio.loadmat('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/data/U_span.mat', squeeze_me=True)
U_span = mat['U_span']

mat = spio.loadmat('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/data/Rij_span.mat', squeeze_me=True)
Rij_span = mat['Rij_span']
uu_span = Rij_span[0:5,:]

mat = spio.loadmat('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/data/Lij_span.mat', squeeze_me=True)
Lij_span = mat['Lij_span']
xLu_span = Lij_span[0:5,:]

fig = '/home/giacomol/Desktop/'
matplotlib.rcParams.update({'font.size': 20})
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True
				   
leg2 = leg
leg2.append('exp')
leg3 = leg2
leg3.append('Inlet')

# BUILDING --------------------------------------------------------------------
f_max = interp1d(zexp, numpy.min(U_span, axis=0), kind='linear')
f_min = interp1d(zexp, numpy.max(U_span, axis=0), kind='linear')

z = numpy.linspace(min(zexp), max(zexp), 10)

plt.figure(figsize=(6,4))
matplotlib.rcParams.update({'font.size': 25})
plt.fill_betweenx(z, f_min(z), f_max(z), color='silver')
plt.plot(numpy.mean(U_span, axis=0), zexp, '-.', color='gray')
plt.plot(U[:,0], zprobes, '--', color='black')
plt.plot(U[:,1:], zprobes)
plt.xlabel(r'$U [m/s]$')
plt.ylim([0,4])
plt.xlim([0,10])
plt.ylabel(r'$y [m]$')
plt.tight_layout()
#plt.legend(['exp','base'])
plt.savefig('/home/giacomol/Desktop/U_step2.png')
plt.show()

plt.figure(figsize=(6,4))
matplotlib.rcParams.update({'font.size': 25})
plt.plot(k_exp, zprobes, '-.', color='gray')
plt.fill_betweenx(zprobes, k_exp_min, k_exp_max, color='silver')
plt.plot(k[:,0], zprobes, '--', color='black')
plt.plot(k[:,1:], zprobes)
plt.xlabel(r'$k [m^2/s^2]$')
plt.ylim([0,4])
plt.xlim([0,1.5])
plt.ylabel(r'$y [m]$')
plt.tight_layout()
plt.locator_params(nbins=3)
plt.savefig(os.path.join(fig, 'k_step2.png'))
plt.show()

plt.figure(figsize=(6,4))
matplotlib.rcParams.update({'font.size': 25})
plt.plot(Tu_exp, zprobes, '-.', color='gray')
plt.fill_betweenx(zprobes, Tu_exp_min, Tu_exp_max, color='silver')
plt.plot(Tu[:,0], zprobes, '--', color='black')
plt.plot(Tu[:,1:], zprobes)
plt.xlabel(r'$T_u [s]$')
plt.ylim([0,4])
plt.xlim([0,0.3])
plt.ylabel(r'$y [m]$')
plt.tight_layout()
plt.locator_params(nbins=3)
plt.savefig(os.path.join(fig, 'Tu_step2.png'))
plt.show()

plt.plot(vv_exp, zprobes, '-.', color='gray')	
plt.fill_betweenx(zprobes, vv_exp_min, vv_exp_max, color='silver')
plt.plot(vv_0, z_vv)
plt.tight_layout()
plt.xlabel(r"$\overline{v'^2} [m^2/s^2]$")
plt.ylim([0,4])
plt.xlim([0,1.3])
plt.ylabel(r'$y [m]$')
plt.savefig(os.path.join(fig, 'vv_0_step2.eps'))
plt.show()

plt.plot(ww_exp, zprobes, '-.', color='gray')	
plt.fill_betweenx(zprobes, ww_exp_min, ww_exp_max, color='silver')
plt.plot(ww_0, z_ww)
plt.tight_layout()
plt.xlabel(r"$\overline{w'^2} [m^2/s^2]$")
plt.ylim([0,4])
plt.xlim([0,1.6])
plt.ylabel(r'$y [m]$')
plt.savefig(os.path.join(fig, 'ww_0_step2.eps'))
plt.show()

plt.plot(uv_exp, zsample, '-.', color='gray')	
plt.fill_betweenx(zsample, uv_exp_min, uv_exp_max, color='silver')
plt.plot(uv, zprobes)
plt.gcf().subplots_adjust(bottom=0.2)
plt.xlabel(r"$\overline{u'v'} [m^2/s^2]$")
plt.ylim([0,4])
plt.xlim([-0.4,0.1])
plt.ylabel(r'$y [m]$')
#plt.savefig(os.path.join(fig, 'uv_b.eps'))
plt.show()

matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize=(4,6))
plt.plot(xLu_exp, zprobes, '-.', color='gray')	
plt.fill_betweenx(zprobes, xLu_exp_min, xLu_exp_max, color='silver')
plt.plot(xLu[:,0], zprobes, '--k', linewidth=2.0)	
plt.plot(xLu[:,1:], zprobes)	
plt.gcf().subplots_adjust(left=0.3)
plt.gcf().subplots_adjust(bottom=0.2)
plt.ylabel(r'$y [m]$')
plt.xlabel(r'$^xL_u [m]$')
plt.ylim([0,4])
plt.xlim([0,2.5])
plt.savefig(os.path.join(fig, 'xLu_sens.eps'))
plt.show()

plt.figure(figsize=(4,6))
plt.plot(xLv_exp, zprobes, '-.', color='gray')	
plt.fill_betweenx(zprobes, xLv_exp_min, xLv_exp_max, color='silver')
plt.plot(xLv[:,0], zprobes, '--k', linewidth=2.0)
plt.plot(xLv[:,1:], zprobes)		
plt.gcf().subplots_adjust(left=0.3)
plt.gcf().subplots_adjust(bottom=0.2)
plt.ylabel(r'$y [m]$')
plt.xlabel(r'$^xL_v [m]$')
plt.ylim([0,4])
plt.xlim([0,0.6])
plt.savefig(os.path.join(fig, 'xLv_sens.eps'))
plt.show()

plt.figure(figsize=(6.5,6))
plt.plot(xLw_exp, zprobes, '-.', color='gray')	
plt.fill_betweenx(zprobes, xLw_exp_min, xLw_exp_max, color='silver')
plt.plot(xLw[:,0], zprobes, '--k', linewidth=2.0)
plt.plot(xLw[:,1:], zprobes)	
plt.gcf().subplots_adjust(left=0.2)
plt.gcf().subplots_adjust(bottom=0.2)
plt.gcf().subplots_adjust(right=0.55)
plt.xlabel(r'$^xL_w [m]$')
plt.ylabel(r'$y [m]$')
plt.xlim([0,1.0])
plt.ylim([0,4])
plt.legend(['exp','base',r"$\overline{u'^2}$",r"$\overline{v'^2}$",
		   r"$\overline{w'^2}$",r"$T_u$",r"$T_v$",r"$T_w$"], 
		   bbox_to_anchor=(1.0, 1.05), ncol=1)
plt.savefig(os.path.join(fig, 'xLw_sens.eps'))
plt.show()


# CELLZONE --------------------------------------------------------------------
plt.plot(vv_exp, zsample, '-.', color='gray')	
plt.fill_betweenx(zsample, vv_exp_min, vv_exp_max, color='silver')
plt.plot(vv_0, z_vv)	
#plt.plot(vv_P[:,0,0], vv_P[:,1,0], 'ok')
plt.legend(['exp', 'step0', 'step1', 'step2'])
plt.gcf().subplots_adjust(bottom=0.2)
plt.xlabel(r"$\overline{v'^2} [m^2/s^2]$")
plt.ylabel(r'$y [m]$')
plt.locator_params(nbins=4)
#plt.savefig(os.path.join(fig, 'vv_0.eps'))
plt.show()

plt.plot(ww_exp, zsample, '-.', color='gray')	
plt.fill_betweenx(zsample, ww_exp_min, ww_exp_max, color='silver')
plt.plot(ww_0, z_ww)	
#plt.plot(ww_P[:,0,0], ww_P[:,1,0], 'ok')
plt.gcf().subplots_adjust(bottom=0.2)
plt.xlabel(r"$\overline{w'^2} [m^2/s^2]$")
plt.ylabel(r'$y [m]$')
plt.locator_params(nbins=4)
#plt.savefig(os.path.join(fig, 'ww_0.eps'))
plt.show()

plt.plot(vv_exp, zsample, '-.', color='gray')	
plt.fill_betweenx(zsample, vv_exp_min, vv_exp_max, color='silver')
plt.plot(vv_0[:,0], z_vv[:,0])
#plt.plot(vv_P[:,0,0], vv_P[:,1,0], 'ok')
plt.gcf().subplots_adjust(bottom=0.2)
plt.xlabel(r"$\overline{v'^2} [m^2/s^2]$")
plt.locator_params(nbins=4)
plt.ylabel(r'$y [m]$')
plt.legend(['exp', 'Bezier curve', 'control points'])
#plt.savefig(os.path.join(fig, 'vv_bezier.eps'))
plt.show()

plt.plot(ww_exp, zsample, '-.', color='gray')	
plt.fill_betweenx(zsample, ww_exp_min, ww_exp_max, color='silver')
plt.plot(ww_0[:,0], z_ww[:,0])
plt.plot(ww_P[:,0,0], ww_P[:,1,0], 'ok')
plt.gcf().subplots_adjust(bottom=0.2)
plt.xlabel(r"$\overline{w'^2} [m^2/s^2]$")
plt.locator_params(nbins=4)
plt.ylabel(r'$y [m]$')
#plt.savefig(os.path.join(fig, 'ww_bezier.eps'))
plt.show()
