import numpy as np
import pandas as pd

import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, Subset
from torchvision import transforms

Nfeatures = 8
Nclasses = 4
height = 8
width = 2

class Dataset(Dataset):
    def __init__(self, path, split, transform, use_gpu):
        super().__init__()
        self.use_gpu = use_gpu
        self.features = []
        self.labels = []
        self.transform = transform

        np.random.seed(42)
        
        if split == 'train' or split == 'valid':
                        
            images_A = np.load(path + '/train/images.npy')
            labels_A = np.load(path + '/train/labels.npy')
                        
            images_D = np.load(path + '/valid/images.npy')
            labels_D = np.load(path + '/valid/labels.npy')
            
            images = np.concatenate((images_A, images_D), axis=0)
            labels = np.concatenate((labels_A, labels_D), axis=0)

            if self.transform:
                
                for i in range(images.shape[0]):
                                        
                    rotated_images = np.zeros_like(images)
                    for j in range(images.shape[2]):
                        transform_image = transforms.ToPILImage()(np.float32(images[i,:,j]))
                        rotated_images[i,:,j] = transforms.functional.vflip(transform_image)
                                    
                images = np.concatenate((images, rotated_images), axis=0)
                labels = np.concatenate((labels, labels), axis=0)
                
            Nim = len(labels)

            index_train = np.random.choice(Nim, size=Nim, replace=False)
            index_valid = index_train[labels >= 2]
            #index_valid = np.setdiff1d(np.arange(Nim), index_train)
            
            if split == 'train':
                index = index_train
            if split == 'valid':
                index = index_valid
            
            self.features = np.transpose(images[index,:,:,images.shape[-1]-Nfeatures:], (0,3,1,2))       
            self.labels = labels[index]
            
        if split == 'test':
                        
            images = np.load(path + '/' + split + '/images.npy')
            labels = np.load(path + '/' + split + '/labels.npy')
    
            self.features = np.transpose(images[:,:,:,images.shape[-1]-Nfeatures:], (0,3,1,2))        
            self.labels = labels
            
        if split == 'tileA':
                        
            images = np.load(path + '/train/images.npy')
            labels = np.load(path + '/train/labels.npy')
    
            self.features = np.transpose(images[:,:,:,images.shape[-1]-Nfeatures:], (0,3,1,2))        
            self.labels = labels
            
        if split == 'tileB':
                        
            images = np.load(path + '/test/images.npy')
            labels = np.load(path + '/test/labels.npy')
    
            self.features = np.transpose(images[:,:,:,images.shape[-1]-Nfeatures:], (0,3,1,2))        
            self.labels = labels            
            
        if split == 'tileD':
                        
            images = np.load(path + '/valid/images.npy')
            labels = np.load(path + '/valid/labels.npy')
    
            self.features = np.transpose(images[:,:,:,images.shape[-1]-Nfeatures:], (0,3,1,2))        
            self.labels = labels
            
    def __getitem__(self, index):

        np.random.seed(42)
                    
        labels = torch.FloatTensor([self.labels[index]])
        features = torch.FloatTensor([self.features[index]])
                
        return features, labels

    def __len__(self):
        return len(self.labels)


def load_data(path, use_gpu=False):

    data_train = Dataset(path, 'train', transform=False, use_gpu=False)
    data_valid = Dataset(path, 'valid', transform=False, use_gpu=False)
    data_test = Dataset(path, 'test', transform=False, use_gpu=False)

    data_A = Dataset(path, 'tileA', transform=False, use_gpu=False)
    data_B = Dataset(path, 'tileB', transform=False, use_gpu=False)
    data_D = Dataset(path, 'tileD', transform=False, use_gpu=False)
        
    return data_train, data_valid, data_test, data_A, data_B, data_D
     