import numpy as np
import pandas as pd
import torch

from utils import plot_panel
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA
import torchvision.transforms as transforms

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
        
    # load coordinates
    coords = np.loadtxt('data/RANS/p_' + angle + 'deg.raw')[:,:3]
    
    m = coords.shape[0]
    n = 8
    
    p_mean = np.loadtxt('data/RANS/p_' + angle + 'deg.raw')[:,3]
    tke = np.loadtxt('data/RANS/k_' + angle + 'deg.raw')[:,3]
    nut = np.loadtxt('data/RANS/nut_' + angle + 'deg.raw')[:,3]
    yPlus = np.loadtxt('data/RANS/yPlus_' + angle + 'deg.raw')[:,3]
    gradp = np.loadtxt('data/RANS/gradp_' + angle + 'deg.raw')[:,3:]
    gradU = np.reshape(np.loadtxt('data/RANS/gradU_' + angle + 'deg.raw')[:,3:], (m,3,3))
    wst = np.loadtxt('data/RANS/wallShearStress_' + angle + 'deg.raw')[:,3]
    
    # compute features ------------------------------------------------------------
    features = np.zeros((m,n))
    
    # 1. mean pressure coefficient
    cp_mean = p_mean/(0.5*Uref**2)
    features[:,0] = cp_mean
            
    # 2. non-dimensional local TKE
    features[:,1] = tke/(Uref**2)
    
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
    features = (features - np.mean(features, axis=0))/np.std(features, axis=0)
    
    # write dataframe -------------------------------------------------------------
    # define data-frame
    faces = np.zeros(m)
    faces[coords[:,0] <= 1e-3] = 1
    faces[coords[:,1] == 2] = 2
    faces[coords[:,2] <= 1e-3] = 5 # switch because of mistake left/right
    faces[coords[:,0] >= 1-1e-3] = 4
    faces[coords[:,2] >= 0.3-1e-3] = 3 # switch because of mistake left/right
    
    dx = 0.025
         
    edges = np.ones(m)
    edges[(faces == 1) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 2) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,0] < X-dx) & (coords[:,0] > dx)] = 0
    edges[(faces == 3) & (coords[:,0] < X-dx) & (coords[:,0] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 4) & (coords[:,2] < Z-dx) & (coords[:,2] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    edges[(faces == 5) & (coords[:,0] < X-dx) & (coords[:,0] > dx) & (coords[:,1] < Y-dx) & (coords[:,1] > dx)] = 0
    
    data = np.column_stack((coords, faces, edges, features))   
    dataFrame = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'face', 'edge',
                                            'cp_mean', 'tke', 'U_inflow',
                                            'yPlus', 'time-scale', 'nut', 'gradp', 'cf'
                                            ])
        
    return features, dataFrame

def data_panel(dataFrame, face, angle, tile):
    
    images = []
    labels = []
    
    all_panels = np.loadtxt('data/windTunnel/designPressure_tile' + tile)
    panels = all_panels[(all_panels[:,4] == angle) & (all_panels[:,5] == face)]
    
    for panel in panels:
        
        if tile == 'A' or tile == 'B':
        
            data = dataFrame[(dataFrame['X'] >= panel[0]) & (dataFrame['X'] <= panel[1]) &
                             (dataFrame['Y'] >= panel[2]) & (dataFrame['Y'] <= panel[3]) &
                             (dataFrame['face'] == face)]
            
            Y_values = np.unique(data['Y'])
            X_values = np.unique(data['X'])
            
            image = np.zeros((len(Y_values), len(X_values), data.values.shape[1]))
            
            for i, Y in enumerate(Y_values):
                image[i,:,:] = data[data['Y'] == Y]
                
        elif tile == 'D':
        
            data = dataFrame[(dataFrame['Z'] >= panel[0]) & (dataFrame['Z'] <= panel[1]) &
                             (dataFrame['Y'] >= panel[2]) & (dataFrame['Y'] <= panel[3]) &
                             (dataFrame['face'] == face)]            
            
            Y_values = np.unique(data['Y'])
            X_values = np.unique(data['Z'])[1:3] # fron/back face have more cells...
            
            image = np.zeros((len(Y_values), len(X_values), data.values.shape[1]))
            
            for i, Y in enumerate(Y_values):
                image[i,0,:] = data[(data['Z'] == X_values[0]) & (data['Y'] == Y)] # fron/back face have more cells...
                image[i,1,:] = data[(data['Z'] == X_values[1]) & (data['Y'] == Y)] # fron/back face have more cells...
                            
        images.append(image)            
        labels.append(int(panel[-1]))
        
    return images, labels

# compute features    
features_00deg, dataFrame_00deg = compute_features('00')
features_10deg, dataFrame_10deg = compute_features('10')
features_20deg, dataFrame_20deg = compute_features('20')
features_30deg, dataFrame_30deg = compute_features('30')
features_45deg, dataFrame_45deg = compute_features('45')

# generate images 
images_A3_00deg, labels_A3_00deg = data_panel(dataFrame_00deg, 3, 0, 'A')
images_B3_00deg, labels_B3_00deg = data_panel(dataFrame_00deg, 3, 0, 'B')
images_D4_00deg, labels_D4_00deg = data_panel(dataFrame_00deg, 4, 0, 'D')
images_D1_00deg, labels_D1_00deg = data_panel(dataFrame_00deg, 1, 0, 'D')

images_A3_10deg, labels_A3_10deg = data_panel(dataFrame_10deg, 3, 10, 'A')
images_B3_10deg, labels_B3_10deg = data_panel(dataFrame_10deg, 3, 10, 'B')
images_D4_10deg, labels_D4_10deg = data_panel(dataFrame_10deg, 4, 10, 'D')

images_A5_10deg, labels_A5_10deg = data_panel(dataFrame_10deg, 5, 10, 'A')
images_B5_10deg, labels_B5_10deg = data_panel(dataFrame_10deg, 5, 10, 'B')
images_D1_10deg, labels_D1_10deg = data_panel(dataFrame_10deg, 1, 10, 'D')

images_A3_20deg, labels_A3_20deg = data_panel(dataFrame_20deg, 3, 20, 'A')
images_B3_20deg, labels_B3_20deg = data_panel(dataFrame_20deg, 3, 20, 'B')
images_D4_20deg, labels_D4_20deg = data_panel(dataFrame_20deg, 4, 20, 'D')

images_A5_20deg, labels_A5_20deg = data_panel(dataFrame_20deg, 5, 20, 'A')
images_B5_20deg, labels_B5_20deg = data_panel(dataFrame_20deg, 5, 20, 'B')
images_D1_20deg, labels_D1_20deg = data_panel(dataFrame_20deg, 1, 20, 'D')

images_A3_30deg, labels_A3_30deg = data_panel(dataFrame_30deg, 3, 30, 'A')
images_B3_30deg, labels_B3_30deg = data_panel(dataFrame_30deg, 3, 30, 'B')
images_D4_30deg, labels_D4_30deg = data_panel(dataFrame_30deg, 4, 30, 'D')

images_A5_30deg, labels_A5_30deg = data_panel(dataFrame_30deg, 5, 30, 'A')
images_B5_30deg, labels_B5_30deg = data_panel(dataFrame_30deg, 5, 30, 'B')
images_D1_30deg, labels_D1_30deg = data_panel(dataFrame_30deg, 1, 30, 'D')

images_A3_45deg, labels_A3_45deg = data_panel(dataFrame_45deg, 3, 45, 'A')
images_B3_45deg, labels_B3_45deg = data_panel(dataFrame_45deg, 3, 45, 'B')
images_D4_45deg, labels_D4_45deg = data_panel(dataFrame_45deg, 4, 45, 'D')

images_A5_45deg, labels_A5_45deg = data_panel(dataFrame_45deg, 5, 45, 'A')
images_B5_45deg, labels_B5_45deg = data_panel(dataFrame_45deg, 5, 45, 'B')
images_D1_45deg, labels_D1_45deg = data_panel(dataFrame_45deg, 1, 45, 'D')

images_train = images_A3_00deg + \
               images_A3_10deg + \
               images_A5_10deg + \
               images_A3_20deg + \
               images_A5_20deg + \
               images_A3_30deg + \
               images_A5_30deg + \
               images_A3_45deg + \
               images_A5_45deg

images_valid = images_D4_00deg + \
               images_D1_00deg + \
               images_D4_10deg + \
               images_D1_10deg + \
               images_D4_20deg + \
               images_D1_20deg + \
               images_D4_30deg + \
               images_D1_30deg + \
               images_D4_45deg + \
               images_D1_45deg      

images_test = images_B3_00deg + \
              images_B3_10deg + \
              images_B5_10deg + \
              images_B3_20deg + \
              images_B5_20deg + \
              images_B3_30deg + \
              images_B5_30deg + \
              images_B3_45deg + \
              images_B5_45deg

labels_train = labels_A3_00deg + \
               labels_A3_10deg + \
               labels_A5_10deg + \
               labels_A3_20deg + \
               labels_A5_20deg + \
               labels_A3_30deg + \
               labels_A5_30deg + \
               labels_A3_45deg + \
               labels_A5_45deg
               
labels_valid = labels_D4_00deg + \
               labels_D1_00deg + \
               labels_D4_10deg + \
               labels_D1_10deg + \
               labels_D4_20deg + \
               labels_D1_20deg + \
               labels_D4_30deg + \
               labels_D1_30deg + \
               labels_D4_45deg + \
               labels_D1_45deg
               
labels_test = labels_B3_00deg + \
              labels_B3_10deg + \
              labels_B5_10deg + \
              labels_B3_20deg + \
              labels_B5_20deg + \
              labels_B3_30deg + \
              labels_B5_30deg + \
              labels_B3_45deg + \
              labels_B5_45deg
'''              
np.save('images_train', images_train)
np.save('labels_train', labels_train)

np.save('images_test', images_test)
np.save('labels_test', labels_test)

np.save('images_valid', images_valid)
np.save('labels_valid', labels_valid)
'''
