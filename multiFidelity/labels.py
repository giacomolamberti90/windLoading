from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from utils import plot_contour
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator 
from scipy.interpolate import interpn
from scipy.interpolate import griddata

#select quadrature point
def load_charLES(workdir):
    # quadrature points
    angles = np.array([ 0.2775616554, 1.777892908, 5.019334521, 10.14314988,  17.0503674,     
                        25.45903128, 34.94759911, 45, 55.05240089, 64.54096872,               
                        72.9496326, 79.85685012, 84.98066548, 88.22210709, 89.72243834])
    
    #Load coords of interest
    coords_RANS = np.loadtxt('data/train/RANS_mesh/p_45deg.raw')[:,:3]
    
    #Load LES-raw data
    back = np.loadtxt('../charLES/workdir.' + str(workdir) + '/back_p_avg.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    front = np.loadtxt('../charLES/workdir.' + str(workdir) + '/front_p_avg.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    left = np.loadtxt('../charLES/workdir.' + str(workdir) + '/left_p_avg.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    right = np.loadtxt('../charLES/workdir.' + str(workdir) + '/right_p_avg.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    top = np.loadtxt('../charLES/workdir.' + str(workdir) + '/top_p_avg.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    
    coords_LES = np.concatenate((back[::1,:3], front[::1,:3], left[::1,:3], right[::1,:3], top[::1,:3]), axis=0)
    p_mean_LES = np.concatenate((back[::1,3], front[::1,3], left[::1,3], right[::1,3], top[::1,3]), axis=0)
    
    back = np.loadtxt('../charLES/workdir.' + str(workdir) + '/back_p_rms.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    front = np.loadtxt('../charLES/workdir.' + str(workdir) + '/front_p_rms.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    left = np.loadtxt('../charLES/workdir.' + str(workdir) + '/left_p_rms.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    right = np.loadtxt('../charLES/workdir.' + str(workdir) + '/right_p_rms.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    top = np.loadtxt('../charLES/workdir.' + str(workdir) + '/top_p_rms.00020000.raw_values.dat', skiprows=2, usecols=(4,5,6,7))
    
    coords_LES = np.concatenate((back[::1,:3], front[::1,:3], left[::1,:3], right[::1,:3], top[::1,:3]), axis=0)
    p_rms_LES = np.concatenate((back[::1,3], front[::1,3], left[::1,3], right[::1,3], top[::1,3]), axis=0)
    
    # Translation of the coords around the origin
    transl1 = np.array([0.5, 0, 0])
    transl2 = np.array([0.5, 0, 0.15])
    
    coords_LES_transl = coords_LES - transl1
    
    # rotate coordinates
    ang = float(angles[workdir-1]) * np.pi/180
    rotation_matrix = np.array([[np.cos(ang), 0, np.sin(ang)],
                                [0, 1, 0],
                                [-np.sin(ang), 0, np.cos(ang)]])
    
    
    coords_LES_rotated = coords_LES_transl.dot(rotation_matrix)
    
    # Translation of coords in same location as RANS
    coords_LES_fix = coords_LES_rotated + transl2
    
    # We need to interpolate to get the data at the RANS points
    # 'Nearest' takes the value closer to the point of interest (fast)
    # 'Linear' takes more time if you use many points
    p_mean_RANS = griddata(coords_LES_fix, p_mean_LES, coords_RANS, method='nearest')
    p_rms_RANS = griddata(coords_LES_fix, p_rms_LES, coords_RANS, method='nearest')
    #np.savetxt('Cp_0.out', Cpi, delimiter=' ') # X is an array
    
    cp_mean_RANS = p_mean_RANS/(0.5*7.7**2)
    cp_rms_RANS = p_rms_RANS/(0.5*7.7**2)
    
    #np.savetxt('../charLES/workdir.' + str(workdir) + '/cp_mean_LES.out', cp_mean_RANS[::5], delimiter='')
    #np.savetxt('../charLES/workdir.' + str(workdir) + '/cp_rms_LES.out', cp_rms_RANS[::5], delimiter='')
    
    #Plot Cp on points of interest and Cp with all points from LES
    #plot_contour(cp_mean_RANS[::5], 2, [-1.5, 1.5], "")
    #plot_contour(cp_rms_RANS[::5], 5, [0, 0.4], "", "angle = " + str(angles[workdir-1])[:3])
    

#select quadrature point
def load_myLES(workdir):
    # quadrature points
    angles = ['10', '45', '80']
        
    #Load LES-raw data
    cp_mean = np.loadtxt('../myLES/pMean_' + angles[workdir-1] + 'deg.raw')/(0.5*7.7**2)
    cp_rms = np.sqrt(np.loadtxt('../myLES/pPrime2Mean_' + angles[workdir-1] + 'deg.raw'))/(0.5*7.7**2)
    
    #cp_rms = np.loadtxt('../RANS/workdir.' + str(workdir) + '/cp_rms_linear.txt')
    
    np.savetxt('../7order/workdir.' + str(workdir) + '/cp_mean_LES.out', cp_mean, delimiter='')
    np.savetxt('../7order/workdir.' + str(workdir) + '/cp_rms_LES.out', cp_rms, delimiter='')
    
    #np.savetxt('../7order/workdir.' + str(workdir) + '/cp_rms_LR_RANS.out', cp_rms[::5], delimiter='')
    
if __name__ == "__main__":
    
    for workdir in range(1,4):
        print(workdir)
        load_myLES(workdir)
    