import numpy as np
import torch

from torch.autograd import Variable
from torch.utils.data import DataLoader

from loader import load_data
from model import CNNNet
from tqdm import tqdm

from utils import plot_class, plot_confusion_matrix
from sklearn.metrics import roc_auc_score

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
        class_weights = torch.tensor([1/548, 1/329, 1/19, 1/4])
        #class_weights = torch.tensor([1/548, 1/329, 1/4, 1/4])

        if train:
            optimizer.zero_grad()

        if loader.dataset.use_gpu:
            features = features.cuda()
            labels = labels.cuda()
            class_weights = class_weights.cuda()
            
        features = Variable(features)
        labels = Variable(labels)
        
        preds = model.forward(features)
                
        loss = torch.nn.CrossEntropyLoss(weight=class_weights)
        output = loss(preds, torch.squeeze(labels, dim=1).long())
                
        total_loss += output.item()
                
        pred_npy = np.argmax(preds.data.cpu().numpy(), axis=1)
        label_npy = labels.data.cpu().numpy()
                
        preds_list.append(np.squeeze(pred_npy))
        labels_list.append(np.squeeze(label_npy))

        #preds.register_hook(lambda grad: print(grad))
        
        if train:
            output.backward()
            optimizer.step()
            
        num_batches += 1
        
    avg_loss = total_loss / num_batches
    
    labels_npy = np.concatenate(labels_list, axis=0)
    preds_npy = np.concatenate(preds_list, axis=0)

    acc = np.mean((preds_npy == labels_npy))
        
    return avg_loss, acc, preds_npy, labels_npy

def evaluate(path, split, angle, face, model_path, use_gpu, filename):
    
    data_train, data_valid, data_test, data_A, data_B, data_D = load_data(path)
        
    model = CNNNet()
    state_dict = torch.load(model_path, map_location=(None if use_gpu else 'cpu'))
    model.load_state_dict(state_dict)

    if use_gpu:
        model = model.cuda()

    if split == 'tileA':        
        loader = DataLoader(data_A, batch_size=32, num_workers=12, shuffle=False)

    elif split == 'tileB':        
        loader = DataLoader(data_B, batch_size=32, num_workers=12, shuffle=False)
        
    elif split == 'tileD':
        loader = DataLoader(data_D, batch_size=32, num_workers=12, shuffle=False)
        
    else:
        raise ValueError("split must be 'train', 'valid', or 'test'")
        
    _, _, preds, labels = run_model(model, loader)
   
    # figure
    acc = plot_class(split, angle, face, preds, labels, filename)
      
    return preds, labels, acc
 
if __name__ == '__main__':
    
    #weights = 'best_models/weighted_ce_loss_wd0p01_valid2more'
    weights = 'models/val0.854167_train0.692308_epoch32'
    
    acc_train = 0; acc_test = 0
    
    # left face at 0deg -------------------------------------------------------
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=0, face=3, model_path=weights, use_gpu=False, filename='design_A3_00')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=0, face=3, model_path=weights, use_gpu=False, filename='design_B3_00')

    acc_train += acc_A; acc_test += acc_B;
    
    # left face at 10deg ------------------------------------------------------
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=10, face=3, model_path=weights, use_gpu=False, filename='design_A3_10')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=10, face=3, model_path=weights, use_gpu=False, filename='design_B3_10')

    acc_train += acc_A; acc_test += acc_B;
    
    # right face at 10deg -----------------------------------------------------
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=10, face=5, model_path=weights, use_gpu=False, filename='design_A5_10')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=10, face=5, model_path=weights, use_gpu=False, filename='design_B5_10')

    acc_train += acc_A; acc_test += acc_B;
    
    # left face at 20deg -------------------------------------------------------  
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=20, face=3, model_path=weights, use_gpu=False, filename='design_A3_20')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=20, face=3, model_path=weights, use_gpu=False, filename='design_B3_20')

    acc_train += acc_A; acc_test += acc_B;
    
    # right face at 20deg -------------------------------------------------------  
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=20, face=5, model_path=weights, use_gpu=False, filename='design_A5_20')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=20, face=5, model_path=weights, use_gpu=False, filename='design_B5_20')

    acc_train += acc_A; acc_test += acc_B;
    
    # left face at 30deg -------------------------------------------------------   
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=30, face=3, model_path=weights, use_gpu=False, filename='design_A3_30')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=30, face=3, model_path=weights, use_gpu=False, filename='design_B3_30')

    acc_train += acc_A; acc_test += acc_B;
    
    # right face at 30deg -------------------------------------------------------    
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=30, face=5, model_path=weights, use_gpu=False, filename='design_A5_30')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=30, face=5, model_path=weights, use_gpu=False, filename='design_B5_30')

    acc_train += acc_A; acc_test += acc_B;
    
    # left face at 45deg -------------------------------------------------------     
    _, _, acc_A = evaluate(path='data/', split='tileA', angle=45, face=3, model_path=weights, use_gpu=False, filename='design_A3_45')
    _, _, acc_B = evaluate(path='data/', split='tileB', angle=45, face=3, model_path=weights, use_gpu=False, filename='design_B3_45')
    
    acc_train += acc_A; acc_test += acc_B;
    
    # right face at 45deg -------------------------------------------------------      
    preds_A, labels_A, acc_A = evaluate(path='data/', split='tileA', angle=45, face=5, model_path=weights, use_gpu=False, filename='design_A5_45')
    preds_B, labels_B, acc_B = evaluate(path='data/', split='tileB', angle=45, face=5, model_path=weights, use_gpu=False, filename='design_B5_45')
    
    acc_train += acc_A; acc_test += acc_B;
    
    acc_train /= 9
    acc_test /= 9
    
    print(f'train accuracy: {acc_train:0.6f}')
    print(f'test accuracy: {acc_test:0.6f}\n')
    
    plot_confusion_matrix(labels_A.astype(int), preds_A.astype(int), np.asarray(['0', '1', '2', '3']), filename='matrix_train', normalize=True)
    plot_confusion_matrix(labels_B.astype(int), preds_B.astype(int), np.asarray(['0', '1', '2', '3']), filename='matrix_test', normalize=True)