import numpy as np

from utils import *
import matplotlib.pyplot as plt

skip = 5
stat = "mean"
clim = [-1.5, 1.5]

# load raw data ===============================================================
coords = np.loadtxt('coords.txt')[::skip,:]
faces = np.loadtxt('faces.txt')[::skip]

cp_RANS = []
cp_LES = []

angles = ['00', '10', '20', '40', '45', '60', '80', '90']

for angle in angles:
    # CharLES
    #cp_RANS.append(np.loadtxt('../RANS/workdir.' + str(workdir) + '/cp_' + stat + '_RANS.out'))
    #cp_LES.append(np.loadtxt('../LES/workdir.' + str(workdir) + '/cp_' + stat + '_LES.out'))
    
    # OpenFOAM
    if stat == "mean":
        cp_RANS.append(np.loadtxt('../RANS/p_' + angle + 'deg.raw')[::skip,-1]/(0.5*7.7**2))
        cp_LES.append(np.loadtxt('../myLES/pMean_' + angle + 'deg.raw')[::skip]/(0.5*7.7**2))
    
    if stat == "rms":
        cp_RANS.append(np.loadtxt('../RANS/cp_rms_PH_' + angle + 'deg.txt')[::skip])
        cp_LES.append(np.sqrt(np.loadtxt('../myLES/pPrime2Mean_' + angle + 'deg.raw')[::skip])/(0.5*7.7**2))
    
    #plot_contour(cp_RANS[workdir-1], skip, clim, 'cpp_mean_contour_LES_3order.png')
    #plot_contour(cp_LES[workdir-1], skip, clim, 'cpp_mean_contour_LES_3order.png')
    
cp_LES = np.asarray(cp_LES)
cp_RANS = np.asarray(cp_RANS)

cp_mean_LES = np.mean(cp_LES, axis=0)
cp_mean_RANS = np.mean(cp_RANS, axis=0)

cp_std_LES = np.std(cp_LES, axis=0)
cp_std_RANS = np.std(cp_RANS, axis=0)

cp_low_LES = np.min(cp_LES, axis=0) #cp_mean_LES - 1.96 * cp_std_LES
cp_low_RANS = np.min(cp_RANS, axis=0) #cp_mean_RANS - 1.96 * cp_std_RANS
  
cp_up_LES = np.max(cp_LES, axis=0) #cp_mean_LES + 1.96 * cp_std_LES
cp_up_RANS = np.max(cp_RANS, axis=0) #cp_mean_RANS + 1.96 * cp_std_RANS

# load data PCE ===============================================================
"""
cp_mean_LES_7order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,0]
cp_low_LES_7order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,1]
cp_up_LES_7order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,2]


cp_mean_LES_5order = np.loadtxt('../5order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,0]
cp_low_LES_5order = np.loadtxt('../5order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,1]
cp_up_LES_5order = np.loadtxt('../5order/skip' + str(skip) + '/PCE_' + stat + '_LES.txt')[:,2]
"""
cp_mean_LES_3order = np.loadtxt('../3order/myLES/skip' + str(skip) + '/PCE_' + stat + '_LES_stdCI.out')[:,0]
cp_low_LES_3order = np.loadtxt('../3order/myLES/skip' + str(skip) + '/PCE_' + stat + '_LES_stdCI.out')[:,1]
cp_up_LES_3order = np.loadtxt('../3order/myLES/skip' + str(skip) + '/PCE_' + stat + '_LES_stdCI.out')[:,2]
cp_angles_LES_3order = np.loadtxt('../3order/myLES/skip' + str(skip) + '/PCE_' + stat + '_LES_angles.txt')

"""
cp_mean_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,0]
cp_low_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat +  '_RANS.out')[:,1]
cp_up_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,2]
"""

if stat == "mean":
    
    cp_mean_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,0]
    cp_low_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,1]
    cp_up_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order.out')[:,2]
    cp_angles_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_multi_3order_angles.txt')
    
    cp_mean_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,0]
    cp_low_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat +  '_RANS.out')[:,1]
    cp_up_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_RANS.out')[:,2]    
    cp_angles_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_RANS_angles.txt')
    
if stat == "rms":
    
    method = 'PH'
    
    cp_mean_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,0]
    cp_low_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,1]
    cp_up_multi = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,2]
    cp_angles_multi= np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_PH_multi_3order_angles.txt')
    
    cp_mean_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_RANS_stdCI.out')[:,0]
    cp_low_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method +  '_RANS_stdCI.out')[:,1]
    cp_up_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_RANS_stdCI.out')[:,2]    
    cp_angles_RANS_7order = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_PH_RANS_angles.txt')
    
    method = 'NN'
    
    cp_mean_multi_ML = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,0]
    cp_low_multi_ML = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,1]
    cp_up_multi_ML = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_stdCI.out')[:,2]
    cp_angles_multi_ML = np.loadtxt('../7order/myLES/skip' + str(skip) + '/PCE_' + stat + '_' + method + '_multi_3order_angles.txt')
    
"""
cp_mean_multi_5order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_multi_5order.txt')[:,0]
cp_low_multi_5order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_multi_5order.txt')[:,1]
cp_up_multi_5order = np.loadtxt('../7order/skip' + str(skip) + '/PCE_' + stat + '_multi_5order.txt')[:,2]
"""

# figures =====================================================================
mse_LES = []
mse_RANS = []
mse_MF = []
mse_MFML = []

for i, angle in enumerate(angles):
        
    if angle not in ['10', '45', '80', '90']:
        
        print(f'angle: {angle:s}')
        
        mse = np.sqrt(np.mean((cp_LES[i] - cp_angles_LES_3order[:,i])**2))
        mse_LES.append(mse)
        print(f'RMSE LES PCE3: {mse:0.4f}')
        
        mse = np.sqrt(np.mean((cp_LES[i] - cp_angles_RANS_7order[:,i])**2))
        mse_RANS.append(mse)
        print(f'RMSE RANS: {mse:0.4f}')    
        
        if stat == "mean":
            mse = np.sqrt(np.mean((cp_LES[i] - cp_angles_multi[:,i])**2))
            mse_MF.append(mse)
            print(f'RMSE MF: {mse:0.4f}\n')
            
        if stat == "rms":
            mse = np.sqrt(np.mean((cp_LES[i] - cp_angles_multi[:,i])**2))
            mse_MF.append(mse)
            print(f'RMSE MF+PH: {mse:0.4f}')
            mse = np.sqrt(np.mean((cp_LES[i] - cp_angles_multi_ML[:,i])**2))
            mse_MFML.append(mse)
            print(f'RMSE MF+ML: {mse:0.4f}\n')
        
        #plot_contourf(cp_LES[i], skip, clim, 'cp_' + stat + '_LES_' + angle + 'deg.png', 'LES')
        #plot_contourf(cp_RANS[i], skip, clim, 'cp_' + stat + '_patersonHolmes_' + angle + 'deg.png', 'RANS')
        #plot_contourf(cp_angles_LES_3order[:,i], skip, clim, 'cp_' + stat + '_mean_contour_LES_3order_' + angle + 'deg.png', 'LES PCE3')
        #plot_contourf(cp_angles_RANS_7order[:,i], skip, clim, 'cp_' + stat + '_mean_contour_RANS_7order_' + angle + 'deg.png', 'RANS PCE7')
        #plot_contourf(cp_angles_multi_ML[:,i], skip, clim, 'cp_' + stat + '_mean_contour_ML_multi_3order_' + angle + 'deg.png', 'MF-ML')
        #plot_contourf(cp_angles_multi[:,i], skip, clim, 'cp_' + stat + '_mean_PH_contour_multi_3order_' + angle + 'deg.png', 'MF-PH')
        
        for plane in ['0', '1', '2', '3']:
            plt.figure()
            plot_profiles(cp_LES[i], cp_LES[i], cp_LES[i], skip, plane, '-.r', clim)
            plot_profiles(cp_RANS[i], cp_RANS[i], cp_RANS[i], skip, plane, '-.b', clim)
            #plot_profiles(cp_angles_LES_3order[:,i], cp_angles_LES_3order[:,i], cp_angles_LES_3order[:,i], skip, plane, '-ok', clim)
            #plot_profiles(cp_angles_RANS_7order[:,i], cp_angles_RANS_7order[:,i], cp_angles_RANS_7order[:,i], skip, plane, '-b', clim)
            plot_profiles(cp_angles_multi[:,i], cp_angles_multi[:,i], cp_angles_multi[:,i], skip, plane, '-k', clim)
            #plot_profiles(cp_angles_multi_ML[:,i], cp_angles_multi_ML[:,i], cp_angles_multi_ML[:,i], skip, plane, '-c', clim)
            #plt.show()
            
            """
            err = np.max(abs(y_LES-y_RANS))
            print(f'RANS max abs error: {err:0.4f}')
            err = np.max(abs(y_LES-y_MF))
            print(f'MF max abs error: {err:0.4f}')
            """
            
            if stat == "mean":
                plt.ylabel(r"$C_P$")
            if stat == "rms":
                plt.ylabel(r"$C_p'$")

            plt.savefig('profile_' + stat + '_' + angle + 'deg_' + plane + '.png')

"""   
print("total:")    
mse = np.mean(mse_LES)
print(f'RMSE LES PCE3: {mse:0.4f}')
mse = np.mean(mse_RANS)
print(f'RMSE RANS: {mse:0.4f}')
mse = np.mean(mse_MF)
print(f'RMSE MF: {mse:0.4f}')
if stat == "rms":
    mse = np.mean(mse_MFML)
    print(f'RMSE MFML: {mse:0.4f}')
"""
#plot_contour(cp_mean_LES, skip, clim, 'cp_' + stat + '_mean_contour_LES.png', 'LES')
#plot_contour(cp_mean_RANS, skip, clim, 'cp_' + stat + '_mean_contour_RANS.png', 'RANS')
#plot_contour(cp_mean_LES_3order, skip, clim, 'cp_' + stat + '_mean_contour_LES_3order.png', 'LES PCE3')
#plot_contour(cp_mean_RANS_7order, skip, clim, 'cp_' + stat + '_mean_contour_RANS_7order.png', 'RANS PCE7')
#plot_contour(cp_mean_multi_ML, skip, clim, 'cp_' + stat + '_mean_contour_ML_multi_3order.png', 'MF-ML')
#plot_contour(cp_mean_multi_PH, skip, clim, 'cp_' + stat + '_mean_contour_PH_multi_3order.png', 'MF-PH')
"""
for plane in ['0', '1', '2', '3']:
    plt.figure()
    #plot_profiles(cp_mean_LES, cp_low_LES, cp_up_LES, skip, plane, '-or', clim)
    #plot_profiles(cp_mean_RANS, cp_low_RANS, cp_up_RANS, skip, plane, '-ok', clim)
    plot_profiles(cp_mean_LES_3order, cp_low_LES_3order, cp_up_LES_3order, skip, plane, '-or', clim)
    plot_profiles(cp_mean_RANS_7order, cp_low_RANS_7order, cp_up_RANS_7order, skip, plane, '-ok', clim)
    plot_profiles(cp_mean_multi, cp_low_multi, cp_up_multi, skip, plane, '-b', clim)
    #plot_profiles(cp_mean_multi_ML, cp_low_multi_ML, cp_up_multi_ML, skip, plane, '-c', clim)
    #plt.show()
    if stat == "mean":
        plt.ylabel(r"$C_P$")
    if stat == "rms":
        plt.ylabel(r"$C_p'$")
    
    if plane == '0':
        plt.legend(['LES: PCE3', 'RANS+PH: PCE7', 'MF'])
    plt.savefig('profile_' + stat + '_' + plane + '.png')
"""