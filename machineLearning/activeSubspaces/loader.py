import numpy as np
import pandas as pd

def load_data(ang):
           
    variables = ['cp_mean', 'tke', 'U_inflow', 'gradp', 'cf']
    
    cpmin = 0.01
    cpmax = 1.0
    
    if ang == 0:
        
        dataFrame = pd.read_csv('../RMSpressure/data/train/RANS_mesh/dataFrame_00deg')
        dataFrame['angle'] = ang * np.ones(len(dataFrame), dtype=np.int32)
                            
    else:
        dataFrame = pd.read_csv('../RMSpressure/data/train/RANS_mesh/dataFrame_' + str(ang) + 'deg')
        dataFrame['angle'] = ang * np.ones(len(dataFrame), dtype=np.int32)

    dataFrame = dataFrame.sample(frac=1, random_state=0).reset_index(drop=True)        
    dataFrame = dataFrame[(dataFrame['cp_rms'] > cpmin) & (dataFrame['cp_rms'] < cpmax)]
    dataFrame = dataFrame[(dataFrame['edge'] == 0)]
    #dataFrame = dataFrame[(dataFrame['face'] == 4)]
    
    labels = dataFrame['cp_rms'].values
    features = dataFrame[variables].values
    #features = (features - np.mean(features, axis=0))/np.std(features, axis=0)
    
    return dataFrame, features, labels