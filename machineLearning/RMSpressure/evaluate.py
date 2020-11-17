import numpy as np
import torch

from torch.autograd import Variable
from torch.utils.data import DataLoader

from loader import load_data_1D, load_data_2D, load_data_tiles, MSELoss
from model import RMSNet, RMSCNNNet
from tqdm import tqdm

from utils import plot_contour, plot_contour_face, plot, plot_profiles
import matplotlib.pyplot as plt
from matplotlib import ticker

def run_model(model, loader, train=False, optimizer=None):
    
    preds_list = []
    labels_list = []

    if train:
        model.train()
    else:
        model.eval()

    total_loss = 0.
    num_batches = 0
        
    
    for batch in loader:
        
        features, labels = batch        
        if train:
            optimizer.zero_grad()

        if loader.dataset.use_gpu:
            features = features.cuda()
            labels = labels.cuda()
            
        features = Variable(features)
        labels = Variable(labels)
        
        preds = model.forward(features)
        
        loss = MSELoss(preds, labels)
        total_loss += loss.item()
                
        pred_npy = preds.data.cpu().numpy()
        label_npy = labels.data.cpu().numpy()
        
        preds_list.append(np.squeeze(pred_npy))
        labels_list.append(np.squeeze(label_npy))

        #preds.register_hook(lambda grad: print(grad))
        
        if train:
            loss.backward()
            optimizer.step()
            
        num_batches += 1
        
    avg_loss = total_loss / num_batches
        
    labels_npy = np.concatenate(labels_list, axis=0)
    preds_npy = np.concatenate(preds_list, axis=0)
    
    mse = np.sqrt(np.mean((preds_npy - labels_npy)**2))
        
    return avg_loss, mse, preds_npy, labels_npy

def evaluate(path, split, model_path, use_gpu):
    
    data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_1D(path)
    #data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_2D(path, 'test')
    #data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_tiles(path)
    
    model = RMSNet()
    state_dict = torch.load(model_path, map_location=(None if use_gpu else 'cpu'))
    model.load_state_dict(state_dict)

    if use_gpu:
        model = model.cuda()

    if split == 'train':
        loader = DataLoader(data_train, batch_size=len(data_train), num_workers=12, shuffle=False)
    elif split == 'valid':
        loader = DataLoader(data_valid, batch_size=len(data_valid), num_workers=12, shuffle=False)        
    elif split == '00deg':
        loader = DataLoader(data_00deg, batch_size=len(data_00deg), num_workers=12, shuffle=False)
    elif split == '20deg':
        loader = DataLoader(data_20deg, batch_size=len(data_20deg), num_workers=12, shuffle=False)
    elif split == '45deg':
        loader = DataLoader(data_45deg, batch_size=len(data_45deg), num_workers=12, shuffle=False)
    else:
        raise ValueError("split must be 'train', 'valid', or 'test'")

    loss, mse, preds, labels = run_model(model, loader)

    print(f'{split} loss: {loss:0.6f}')
    print(f'{split} MSE: {mse:0.6f}\n')
        
    plot(labels, preds, '.r')
    
    if split == '45deg':
        plot_contour(loader.dataset.coords, labels, [0, 0.3], '')
        plot_contour(loader.dataset.coords, loader.dataset.empiricalModels[:,0], [0, 0.3], '')
        plot_contour(loader.dataset.coords, preds, [0, 0.3], '')
        
        plot_profiles(loader.dataset, preds, preds/1000, 'Ynorm')
           
        #plot_contour_face(loader.dataset.coords, labels, [0, 0.3], face, '')
        #plot_contour_face(loader.dataset.coords, loader.dataset.empiricalModels[:,0], [0, 0.3], face, '')
        #plot_contour_face(loader.dataset.coords, preds, [0, 0.3], face, '')
    
    '''
    coords_1D = []
    labels_1D = []
    preds_1D = []
    for i in range(labels.shape[0]):
        coords_1D.append(np.reshape(loader.dataset.coords[i], (48*16,2)))
        labels_1D.append(np.reshape(labels[i], (48*16,1)))
        preds_1D.append(np.reshape(preds[i], (48*16,1)))
    
    coords_1D = np.reshape(np.asarray(coords_1D), (labels.shape[0]*48*16,2))
    labels_1D = np.reshape(np.asarray(labels_1D), (labels.shape[0]*48*16,1))   
    preds_1D = np.reshape(np.asarray(preds_1D), (labels.shape[0]*48*16,1))   
    faces = 5*np.ones_like(preds_1D)
    
    #plot_profiles(coords_1D, labels_1D, preds_1D, faces, 'Ynorm')
    
    # labels
    plt.figure(figsize=(4,5))
    plt.rcParams.update({'font.size': 18})
    for i in range(labels.shape[0]):
        coords = np.reshape(loader.dataset.coords[i], (48*16,2))
        label = np.reshape(labels[i], (48*16,1))
        plt.scatter(coords[:,0], coords[:,1], c=label, marker='.', cmap='hot')
        plt.clim([0, 0.3])
        plt.xlabel(r"$x[m]$"); plt.ylabel(r"$y[m]$")

    cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    plt.xticks([]); plt.yticks([])
    plt.tight_layout()    
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    
    # prediction
    plt.figure(figsize=(4,5))
    plt.rcParams.update({'font.size': 18})
    for i in range(labels.shape[0]):
        coords = np.reshape(loader.dataset.coords[i], (48*16,2))
        pred = np.reshape(preds[i], (48*16,1))
        plt.scatter(coords[:,0], coords[:,1], c=pred, marker='.', cmap='hot')
        plt.clim([0, 0.3])
        plt.xlabel(r"$x[m]$"); plt.ylabel(r"$y[m]$")

    cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    plt.xticks([]); plt.yticks([])
    plt.tight_layout()    
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    '''
    return preds, labels
 
if __name__ == '__main__':
    
    weights = 'models/val0.037031_train0.032620_epoch24'
    evaluate(path='data/', split='train', model_path=weights, use_gpu=False)
    evaluate(path='data/', split='00deg', model_path=weights, use_gpu=False)
    evaluate(path='data/', split='20deg', model_path=weights, use_gpu=False)    
    evaluate(path='data/', split='45deg', model_path=weights, use_gpu=False)
    
    