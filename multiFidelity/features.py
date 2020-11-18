import numpy as np
import pandas as pd
from utils import *

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
    
def compute_features(workdir):
        
    # load coordinates
    coords = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/p_highRiseSampling.raw')[:,:3]
    
    # features
    m = coords.shape[0]
    n = 8
    
    p_mean = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/p_highRiseSampling.raw')[:,3]
    tke = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/k_highRiseSampling.raw')[:,3]
    epsilon = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/epsilon_highRiseSampling.raw')[:,3]
    nut = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/nut_highRiseSampling.raw')[:,3]
    yPlus = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/yPlus_highRiseSampling.raw')[:,3]
    gradp = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/grad(p)_highRiseSampling.raw')[:,3:]
    gradU = np.reshape(np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/grad(U)_highRiseSampling.raw')[:,3:], (m,3,3))
    wst = np.loadtxt('../RANS/workdir.' + str(workdir) + '/postProcessing/surfaces/3000/wallShearStress_highRiseSampling.raw')[:,3]
    
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
    features[:,7] = wst / (0.5*Uref**2)
    
    # normalize features
    #features = (features - np.mean(features, axis=0))/np.std(features, axis=0)
    
    # Paterson Holmes
    cp_rms_PH = (tke/3 + 0.816*abs(cp_mean)*U*np.sqrt(kref))/(0.5*Uref**2)

    # output
    cp_rms = np.ones((m,1))#np.loadtxt('../LES/workdir.' + str(2) + '/cp_rms.out')
    
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
    
    #np.savetxt('coords.txt', coords, delimiter=' ')
    #np.savetxt('edges.txt', edges, delimiter=' ')
    
    np.savetxt('../7order/workdir.' + str(workdir) + '/cp_mean_RANS.out', cp_mean, delimiter='')
    np.savetxt('../7order/workdir.' + str(workdir) + '/cp_rms_PH_RANS.out', cp_rms_PH, delimiter='')
    #dataFrame.to_csv('../RANS/workdir.' + str(workdir) + '/dataFrame')
        
    return features, dataFrame

for workdir in range(1,2):
    print(workdir)
    features, dataFrame = compute_features(workdir)