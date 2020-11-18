import numpy as np

from utils import *
import matplotlib.pyplot as plt

skip = 5
stat = "rms"
clim = [0, 0.5]

# load raw data ===============================================================
coords = np.loadtxt('coords.txt')[::skip,:]
faces = np.loadtxt('faces.txt')[::skip]

cp_RANS = []
cp_LES = []

# OpenFOAM
#angles = ['00', '10', '20', '40', '45', '60', '80', '90']

# CharLES
angles = np.array([ 0.2775616554, 1.777892908, 5.019334521, 10.14314988,  17.0503674,     
                    25.45903128, 34.94759911, 45, 55.05240089, 64.54096872,               
                    72.9496326, 79.85685012, 84.98066548, 88.22210709, 89.72243834])

# load data PCE ===============================================================
cp_mean_LES_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,0]
cp_low_LES_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,1]
cp_up_LES_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,2]

cp_mean_LES_3order = np.loadtxt('../3order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,0]
cp_low_LES_3order = np.loadtxt('../3order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,1]
cp_up_LES_3order = np.loadtxt('../3order/charLES/skip' + str(skip) + '/PCE_' + stat + '_LES.out')[:,2]

"""
cp_mean_RANS_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,0]
cp_low_RANS_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat +  '_RANS.out')[:,1]
cp_up_RANS_7order = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,2]
"""
    
cp_mean_multi = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,0]
cp_low_multi = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,1]
cp_up_multi = np.loadtxt('../7order/charLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,2]

# figures =====================================================================
plot_contour(cp_mean_LES_3order, skip, clim, 'cp_' + stat + '_mean_contour_LES_3order.png', 'LES PCE3')
plot_contour(cp_mean_LES_7order, skip, clim, 'cp_' + stat + '_mean_contour_RANS_7order.png', 'LES PCE7')
plot_contour(cp_mean_multi, skip, clim, 'cp_' + stat + '_mean_contour_multi_3order.png', 'MF')

for plane in ['0', '1', '2', '3']:
    plt.figure()
    plot_profiles(cp_mean_LES_3order, cp_low_LES_3order, cp_up_LES_3order, skip, plane, '-or', clim)
    plot_profiles(cp_mean_LES_7order, cp_low_LES_7order, cp_up_LES_7order, skip, plane, '-ok', clim)
    plot_profiles(cp_mean_multi, cp_low_multi, cp_up_multi, skip, plane, '-b', clim)
    plt.show()
    if stat == "mean":
        plt.ylabel(r"$C_P$")
    if stat == "rms":
        plt.ylabel(r"$C_p'$")
    
    if plane == '0':
        plt.legend(['LES3', 'LES7', 'MF'])
    plt.savefig('profile_' + stat + '_' + plane + '.png')
