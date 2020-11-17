import numpy as np
import pandas as pd

from utils import plot_contour
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator

import random

np.random.seed(1)

# load files ------------------------------------------------------------------
angle = '45'
split = 'train'

# load coordinates
coords_RANS = np.loadtxt('data/train/RANS_mesh/p_45deg.raw')[:,:3]
#coords_LES = np.loadtxt('data/' + split + '/LES_mesh/pPrime2Mean_' + angle + 'deg.raw')[:,:3]

#coords_A = np.loadtxt('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/coords_A0')
#coords_B = np.loadtxt('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/coords_B0')

# labels
#cp_rms_LES = np.sqrt(np.loadtxt('data/' + split + '/RANS_mesh/pPrime2Mean_' + '45' + 'deg.raw'))/(0.5*7.7**2)
#cp_rms_RANS = np.sqrt(np.loadtxt('data/' + split + '/RANS_mesh/pPrime2Mean_' + angle + 'deg.raw'))/(0.5*7.7**2)

#coords_LES = np.loadtxt('../../LES/highRise/20deg/coarse/dynamicK/workdir.14/postProcessing/surfaces/10000/p_highRiseSampling.raw')[:,:3]

# rotate coordinates
ang = -int(angle) * np.pi/180
rotation_matrix = np.array([[np.cos(ang), 0, np.sin(ang)],
                            [0, 1, 0],
                            [-np.sin(ang), 0, np.cos(ang)]])

coords_RANS_rotated = coords_RANS.dot(rotation_matrix)

#coords_LES_rotated = np.around(coords_LES.dot(rotation_matrix), decimals=4)

#plot_contour(coords_RANS, cp_rms_LES, [0, 0.3], '')
#plot_contour(coords_RANS, cp_rms_RANS, [0, 0.3], '')
'''
# select probes:
index = (coords_LES_rotated[:,2] > 0.29) * (coords_LES_rotated[:,1] > 1)
coords_LES_rotated = coords_LES_rotated[index]
coords_LES = coords_LES[index]

#index = np.random.choice(np.nonzero(index)[0], size=1000)

Y_values = np.unique(coords_LES_rotated)
X_values = np.unique(coords_LES_rotated)

Y_values = np.around(Y_values[0::10], decimals=4)
X_values = np.around(X_values, decimals=4)

probes = []
for idx, probe in enumerate(coords_LES):
    
    if (coords_LES_rotated[idx, 0] in X_values) and (coords_LES_rotated[idx, 1] in Y_values):
        probes.append(list(np.around(probe, decimals=4)))
        
# make contour plot
plt.rcParams.update({'font.size': 18})   
x = coords_RANS_rotated[index,0]+0.5
y = coords_RANS_rotated[index,1]
z = cp_rms_LES[index]
plt.scatter(x, y, c=z, marker='.', cmap='hot')
plt.clim([0, 0.3])


coords_LES = np.loadtxt('/home/giacomol/Desktop/Research/windLoading/LES/highRise/90deg/postProcessing/surfaces/36.3/pPrime2Mean_highRiseSampling.raw')[0::100,:3]
'''
f = open('coords.txt', 'w')
for i in range(len(coords_RANS_rotated)):
    f.write('%6.4f %6.4f %6.4f\n' % tuple(coords_RANS_rotated[i]))
f.close()
