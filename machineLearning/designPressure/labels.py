import numpy as np
import pandas as pd

from utils import plot_contour
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator

import random

np.random.seed(1)

# load files ------------------------------------------------------------------
angle = '30'
split = 'train'

# load coordinates
coords_RANS = np.loadtxt('data/train/RANS_mesh/p_45deg.raw')[:,:3]
coords_LES = np.loadtxt('data/' + split + '/LES_mesh/pPrime2Mean_' + angle + 'deg.raw')[:,:3]

coords_A = np.loadtxt('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/coords_A0')
coords_B = np.loadtxt('/home/giacomol/Desktop/Research/windLoading/windTunnel/PoliMi/coords_B0')
    
# labels
#cp_rms_LES = np.sqrt(np.loadtxt('data/' + split + '/LES_mesh/pPrime2Mean_' + angle + 'deg.raw')[:,3])/(0.5*7.7**2)
#cp_rms_RANS = np.sqrt(np.loadtxt('data/' + split + '/RANS_mesh/pPrime2Mean_' + angle + 'deg.raw'))/(0.5*7.7**2)

# rotate coordinates
ang = -int(angle)*np.pi/180
rotation_matrix = np.array([[np.cos(ang), 0, np.sin(ang)],
                            [0, 1, 0],
                            [-np.sin(ang), 0, np.cos(ang)]])

coords_RANS_rotated = coords_RANS.dot(rotation_matrix)
'''
#plot_contour(coords_LES_rotated, cp_rms_LES, [0, 0.3], '')
plot_contour(coords_RANS, cp_rms_RANS, [0, 0.3], '')

# select probes:
index = (coords_LES_rotated[:,2] > 0.299) * (coords_LES_rotated[:,0] > 0.75) * (coords_LES_rotated[:,1] > 1.75)
#index = np.random.choice(np.nonzero(index)[0], size=1000)

Y_values = np.unique(coords_LES_rotated[index,1])
Y_values = Y_values[0::10]

# make contour plot
plt.rcParams.update({'font.size': 18})   
x = coords_LES_rotated[index,0]+0.5
y = coords_LES_rotated[index,1]
z = cp_rms_LES[index]
plt.scatter(x, y, c=z, marker='.', cmap='hot')
plt.clim([0, 0.3])
'''
f = open('coords.txt', 'w')
for i in range(coords_RANS_rotated.shape[0]):
    f.write('%6.7f %6.7f %6.7f\n' % tuple(coords_RANS_rotated[i]))
f.close()    
