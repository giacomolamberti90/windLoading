import numpy as np
import pandas as pd

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Subset
from torchvision import transforms

from utils import plot_contour

face = 5

class Dataset_1D(Dataset):
    def __init__(self, path, split, use_gpu):
        super().__init__()
        self.use_gpu = False
        self.features = []
        self.labels = []
        self.coords = []
        self.faces = []
                
        variables = ['cp_mean', 'tke', 'U_inflow', 'gradp', 'cf']
        emp_models = ['paterson-homes']
        coordinates = ['X', 'Y', 'Z']
        
        cpmin = 0.01
        cpmax = 1.00
            
        if split == 'train':

            data_10 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_10deg')
            data_30 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_30deg')
            data_50 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_50deg')
            data_70 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_70deg')
            data_90 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_90deg')
            
            dataFrame = data_90
            dataFrame = dataFrame[(dataFrame['cp_rms'] > cpmin) & (dataFrame['cp_rms'] < cpmax)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            #dataFrame = dataFrame[(dataFrame['face'] == 5) | (dataFrame['face'] == 3)]
            #dataFrame = dataFrame.drop(['selvam', 'paterson', 'richards'], axis=1)
            
            data_train = dataFrame.sample(frac=1.00, random_state=0)
                        
            self.coords = data_train[coordinates].values
            #self.labels = np.log(dataFrame.values[:,6])
            self.labels = data_train['cp_rms'].values
            self.empiricalModels = data_train[emp_models].values
            self.features = data_train[variables].values
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = data_train['face']
            
        elif split == 'valid':            
 
            data_45 = pd.read_csv(path + '/train/RANS_mesh/dataFrame_45deg')           
            
            dataFrame = data_45
            dataFrame = dataFrame[(dataFrame['cp_rms'] > cpmin) & (dataFrame['cp_rms'] < cpmax)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            #dataFrame = dataFrame[(dataFrame['face'] == face)]
            #dataFrame = dataFrame.drop(['selvam', 'paterson', 'richards'], axis=1)
                        
            self.coords = dataFrame[coordinates].values
            #self.labels = np.log(dataFrame.values[:,6])
            self.labels = dataFrame['cp_rms'].values
            self.empiricalModels = dataFrame[emp_models].values
            self.features = dataFrame[variables].values
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
            
        else:
                        
            dataFrame = pd.read_csv(path + '/train/RANS_mesh/dataFrame_' + split)
            dataFrame = dataFrame[(dataFrame['cp_rms'] > cpmin) & (dataFrame['cp_rms'] < cpmax)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            #dataFrame = dataFrame[(dataFrame['face'] == face)]         
            #dataFrame = dataFrame.drop(['selvam', 'paterson', 'richards'], axis=1)
            
            self.coords = dataFrame[coordinates].values
            #self.labels = np.log(dataFrame.values[:,6])
            self.labels = dataFrame['cp_rms'].values
            self.empiricalModels = dataFrame[emp_models].values
            self.features = dataFrame[variables].values
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
                   
            
    def __getitem__(self, index):
        
        labels = torch.FloatTensor([self.labels[index]])
        features = torch.FloatTensor([self.features[index]])
                
        return features, labels

    def __len__(self):
        return len(self.labels)

class Dataset_2D(Dataset):
    def __init__(self, path, split, transform, use_gpu):
        super().__init__()
        self.use_gpu = use_gpu
        self.features = []
        self.labels = []
        self.coords = []
        self.transform = transform
        
        np.random.seed(42)
        
        if split == '00deg':
            
            images = np.load(path + '/crops/right/images_test_00deg.npy')
            
            self.coords = images[:,:,:,:2]
            self.labels = images[:,:,:,2]
            self.features = images[:,:,:,3:]
            
        if split == '20deg':
            
            images = np.load(path + '/crops/right/images_test_20deg.npy')

            self.coords = images[:,:,:,:2]
            self.labels = images[:,:,:,2]
            self.features = images[:,:,:,3:]
            
        if split == '45deg':
            
            images = np.load(path + '/crops/left/images_45deg.npy')
            
            self.coords = images[:,:,:,:2]
            self.labels = images[:,:,:,2]
            self.features = images[:,:,:,3:]
            
        if split == 'train' or split == 'valid':
            
            images = np.concatenate((np.load(path + '/crops/right/images_20deg.npy'),
                                     np.load(path + '/crops/left/images_20deg.npy')), axis=0)
            
            if self.transform:
                
                for i in range(images.shape[0]):
                    
                    angle = np.random.choice([-20, -10, 10, 20])
                    
                    rotated_images = np.zeros_like(images)
                    for j in range(images.shape[2]):
                        transform_image = transforms.ToPILImage()(np.float32(images[i,:,j]))
                        rotated_images[i,:,j] = transforms.functional.rotate(transform_image, angle)
                                    
                images = np.concatenate((images, rotated_images), axis=0)            
            
            index_train = np.random.choice(len(images), size=len(images)-100, replace=False)
            index_valid = np.setdiff1d(np.arange(len(images)), index_train)
            
            if split == 'train':
                index = index_train
            if split == 'valid':
                index = index_valid
                
            self.coords = images[index,:,:,:2]
            self.labels = images[index,:,:,2]
            self.features = images[index,:,:,3:]
        
    def __getitem__(self, index):

        np.random.seed(42)
                    
        labels = torch.FloatTensor([self.labels[index]])
        features = torch.FloatTensor([self.features[index]])
                
        return features, labels

    def __len__(self):
        return len(self.labels)

class tiles(Dataset):
    def __init__(self, path, split, use_gpu):
        super().__init__()
        self.use_gpu = False
        self.features = []
        self.labels = []
        self.coords = []
        self.faces = []
        
        fin = 5
        fend = 8

        if split == '00deg':
            
            dataFrame = pd.read_csv(path + '/train/RANS_mesh/dataFrame_00deg')
            dataFrame = dataFrame[(dataFrame['cp_rms'] > 0.01)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            
            self.coords = dataFrame.values[:,1:4]
            self.labels = np.log(dataFrame.values[:,6])
            self.empiricalModels = dataFrame.values[:,7:11]
            self.features = dataFrame.values[:,fin+6:fend+6]
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
            
        if split == '20deg':
            
            dataFrame = pd.read_csv(path + '/train/RANS_mesh/dataFrame_20deg') 
            dataFrame = dataFrame[(dataFrame['cp_rms'] > 0.01)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            
            self.coords = dataFrame.values[:,1:4]
            self.labels = np.log(dataFrame.values[:,6])
            self.empiricalModels = dataFrame.values[:,7:11]
            self.features = dataFrame.values[:,fin+6:fend+6]
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
            
        if split == '45deg':
            
            dataFrame = pd.read_csv(path + '/train/RANS_mesh/dataFrame_45deg') 
            dataFrame = dataFrame[(dataFrame['cp_rms'] > 0.01)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            
            self.coords = dataFrame.values[:,1:4]
            self.labels = np.log(dataFrame.values[:,6])
            self.empiricalModels = dataFrame.values[:,7:11]
            self.features = dataFrame.values[:,fin+6:fend+6]
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
            
        if split == 'train':
            
            data_10deg = pd.read_csv(path + '/train/RANS_mesh/dataFrame_10deg')
            data_50deg = pd.read_csv(path + '/train/RANS_mesh/dataFrame_50deg')
            
            dataFrame = data_50deg#.append(data_30deg)
            dataFrame = dataFrame[(dataFrame['rmsCp'] > 0.01)]

            
            data_train = dataFrame #dataFrame.sample(frac=0.9, random_state=0)         
            #data_valid = dataFrame.drop(data_train.index, axis=0)
                        
            self.coords = data_train.values[:,1:3]
            self.labels = np.log(data_train.values[:,4])
            self.features = data_train.values[:,fin:fend]
            
        if split == 'valid':
            
            dataFrame = pd.read_csv(path + '/train/RANS_mesh/dataFrame_45deg')
            dataFrame = dataFrame[(dataFrame['cp_rms'] > 0.01)]
            dataFrame = dataFrame[(dataFrame['edge'] == 0)]
            
            self.coords = dataFrame.values[:,1:4]
            self.labels = np.log(dataFrame.values[:,6])
            self.empiricalModels = dataFrame.values[:,7:11]
            self.features = dataFrame.values[:,fin+6:fend+6]
            self.features = (self.features - np.mean(self.features, axis=0))/np.std(self.features, axis=0)            
            self.faces = dataFrame['face']
            
    def __getitem__(self, index):
        
        labels = torch.FloatTensor([self.labels[index]])
        features = torch.FloatTensor([self.features[index]])
                
        return features, labels

    def __len__(self):
        return len(self.labels)
    
def MSELoss(yhat,y):
    
    loss_MSE = torch.sqrt(torch.mean((yhat-y)**2))
    loss_MAE = torch.mean(torch.abs(yhat-y))
    
    loss = torch.mul(loss_MSE, loss_MSE < 1) + torch.mul(loss_MAE, loss_MAE >= 1)
    
    return loss_MSE

def load_data_1D(path, use_gpu=False):
    
    data_train = Dataset_1D(path, 'train', use_gpu)
    data_valid = Dataset_1D(path, 'valid', use_gpu)
    
    data_10deg = Dataset_1D(path, '10deg', use_gpu)
    data_45deg = Dataset_1D(path, '45deg', use_gpu)
    
    data_00deg = Dataset_1D(path, '00deg', use_gpu)
    data_20deg = Dataset_1D(path, '20deg', use_gpu)
    data_40deg = Dataset_1D(path, '40deg', use_gpu)
    data_60deg = Dataset_1D(path, '60deg', use_gpu)    
    data_80deg = Dataset_1D(path, '80deg', use_gpu)
    
    return data_train, data_valid, data_00deg, data_20deg, data_40deg, data_60deg, data_80deg

def load_data_2D(path, split, use_gpu=False):

    if split == 'test':
        data_train = []
        data_valid = []
        data_00deg = Dataset_2D(path, '00deg', False, use_gpu)
        data_20deg = Dataset_2D(path, '20deg', False, use_gpu)
        data_45deg = Dataset_2D(path, '45deg', False, use_gpu)
    else:
        data_train = Dataset_2D(path, 'train', False, use_gpu)
        data_valid = Dataset_2D(path, 'valid', False, use_gpu)
        data_00deg = []
        data_20deg = []
        data_45deg = []
        
    return data_train, data_valid, data_00deg, data_20deg, data_45deg

def load_data_tiles(path, use_gpu=False):

    data_00deg = tiles(path, '00deg', use_gpu)
    data_20deg = tiles(path, '20deg', use_gpu)
    data_45deg = tiles(path, '45deg', use_gpu)
    
    data_train = tiles(path, 'train', use_gpu)
    data_valid = tiles(path, 'valid', use_gpu)
    
    return data_train, data_valid, data_00deg, data_20deg, data_45deg

