import numpy as np
import pandas as pd
from utils import *
from scipy.interpolate import griddata

# data ------------------------------------------------------------------------
X = 1
Y = 2
Z = 0.3

H = 2 # heigh of building
Uref = 7.7 # velocity at building height
z0 = 0.0032 # roughness length
ustar = Uref*0.4/np.log((H + z0)/z0)
kref = ustar**2/np.sqrt(0.09)
nu = 1e-5
    
def compute_features(angle):
    
    # output
    cp_rms = np.sqrt(np.loadtxt('data/train/RANS_mesh/pPrime2Mean_' + angle + 'deg.raw'))/(0.5*Uref**2)
    
    # load coordinates
    coords = np.loadtxt('data/train/RANS_mesh/p_' + angle + 'deg.raw')[:,:3]
    
    # features
    m = coords.shape[0]
    n = 8
    
    p_mean = np.loadtxt('data/train/RANS_mesh/p_' + angle + 'deg.raw')[:,3]
    tke = np.loadtxt('data/train/RANS_mesh/k_' + angle + 'deg.raw')[:,3]
    epsilon = np.loadtxt('data/train/RANS_mesh/epsilon_' + angle + 'deg.raw')[:,3]
    nut = np.loadtxt('data/train/RANS_mesh/nut_' + angle + 'deg.raw')[:,3]
    yPlus = np.loadtxt('data/train/RANS_mesh/yPlus_' + angle + 'deg.raw')[:,3]
    gradU = np.reshape(np.loadtxt('data/train/RANS_mesh/gradU_' + angle + 'deg.raw')[:,3:], (m,3,3))

    # interpolate
    wst_point = np.loadtxt('data/train/RANS_mesh/wallShearStress_point_' + angle + 'deg.raw')[:,3:]
    
    tauw_point = np.linalg.norm(wst_point, axis=1)
    
    gradp_point = np.loadtxt('data/train/RANS_mesh/gradp_point_' + angle + 'deg.raw')[:,3:]
    coords_point = np.loadtxt('data/train/RANS_mesh/wallShearStress_point_' + angle + 'deg.raw')[:,:3]
    
    tauw = griddata(coords_point, tauw_point, coords, method='nearest')
    gradp = griddata(coords_point, gradp_point, coords, method='nearest')
        
    # compute features ------------------------------------------------------------
    features = np.zeros((m,n))
    
    # 1. mean pressure coefficient
    cp_mean = p_mean/(0.5*Uref**2)
    features[:,0] = cp_mean
            
    # 2. non-dimensional local TKE
    features[:,1] = tke/Uref**2
    
    # 3. non-dimensional velocity at inlet
    U = ustar/0.4 * np.log((coords[:,1] + z0)/z0)
    features[:,2] =  U/Uref
    
    # 4. yPlus 
    features[:,3] = yPlus
    
    # 5. Ratio of turbulence time-scale to mean strain time-scale
    strain = 0.5 * (gradU + np.transpose(gradU, (0,2,1)))
    rotation = 0.5 * (gradU - np.transpose(gradU, (0,2,1)))
    
    rotation_frob2 = np.square(np.linalg.norm(rotation, axis=(1,2)))
    strain_frob2 = np.square(np.linalg.norm(strain, axis=(1,2)))
    
    features[:,4] = strain_frob2*tke/epsilon
    
    # 6. Viscosity ratio
    features[:,5] = nut/(100*nu)
    
    # 7. Ratio of pressure normal stresses to normal shear stresses
    features[:,6] = np.linalg.norm(gradp, axis=1) * H/(0.5*Uref**2)
    
    # 8. wall shear stress
    features[:,7] = tauw / (0.5*Uref**2)
    
    # normalize features
    #features = (features - np.mean(features, axis=0))/np.std(features, axis=0)
    """
    plot_contourf(coords, features[:,0], [np.min(features[:,0]), np.max(features[:,0])], 'cp_mean_00deg.png')
    plot_contourf(coords, features[:,1], [np.min(features[:,1]), np.max(features[:,1])], 'tke_00deg.png')
    plot_contourf(coords, features[:,2], [np.min(features[:,2]), np.max(features[:,2])], 'uinflow_00deg.png')
    plot_contourf(coords, features[:,6], [0, 50], 'gradp_00deg.png')
    plot_contourf(coords, features[:,7], [0, 0.02], 'cf_00deg.png')
    """
    # Paterson Holmes
    cp_rms_PH = (tke/3 + 0.816*abs(cp_mean)*U*np.sqrt(kref))/(0.5*Uref**2)
    
    # write dataframe -------------------------------------------------------------
    # define data-frame
    faces = np.zeros(m)
    faces[coords[:,0] <= 1e-3] = 1
    faces[coords[:,1] == 2] = 2
    faces[coords[:,2] <= 1e-3] = 3
    faces[coords[:,0] >= 1-1e-3] = 4
    faces[coords[:,2] >= 0.3-1e-3] = 5
    
    dx = 0.025
         
    edges = np.ones(m)
    edges[(faces == 1) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 2) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,0] < X-dx) & (coords[:,0] > dx)] = 0
    edges[(faces == 3) & (coords[:,0] < X-dx) & (coords[:,0] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 4) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 5) & (coords[:,0] < X-dx) & (coords[:,0] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    
    data = np.column_stack((coords, faces, edges, cp_rms, features, cp_rms_PH))   
    dataFrame = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'face', 'edge',
                                            'cp_rms', 'cp_mean', 'tke', 'U_inflow',
                                            'yPlus', 'time-scale', 'nut', 'gradp', 'cf',
                                            'paterson-homes'])
    
    dataFrame.to_csv('dataFrame_cf_' + angle + 'deg')
        
    return features, dataFrame

# compute features
for i in range(10):
    features_00deg_new, dataFrame_00deg_new = compute_features('45')
