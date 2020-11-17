import numpy as np
import torch

from datetime import datetime
from pathlib import Path

from evaluate import run_model
from loader import load_data
from model import CNNNet

from torch.utils.data import DataLoader

def train(rundir, path, epochs, learning_rate, use_gpu):
    
    data_train, data_valid, data_test, data_A, data_B, data_D = load_data(path)

    train_loader = DataLoader(data_train, batch_size=32, num_workers=12, shuffle=True)
    valid_loader = DataLoader(data_valid, batch_size=32, num_workers=12, shuffle=True)
    
    model = CNNNet()
    
    if use_gpu:
        model = model.cuda()

    optimizer = torch.optim.Adam(model.parameters(), learning_rate, weight_decay=0.01)
    #scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=.3, threshold=1e-5)

    best_val_acc = -float('inf')

    start_time = datetime.now()
    
    for epoch in range(epochs):
                
        change = datetime.now() - start_time
        print('starting epoch {}. time passed: {}\n'.format(epoch+1, str(change)))
        
        train_loss, train_acc, _, _ = run_model(model, train_loader, train=True, optimizer=optimizer)
        print(f'train loss: {train_loss:0.6f}')
        print(f'train accuracy: {train_acc:0.6f}\n')         
        
        val_loss, val_acc, _, _ = run_model(model, valid_loader)
        print(f'valid loss: {val_loss:0.6f}')
        print(f'valid accuracy: {val_acc:0.6f}\n')
        
        #scheduler.step(val_loss)
        
        if val_acc >= best_val_acc:
            best_val_acc = val_acc

            file_name = f'val{val_acc:0.6f}_train{train_acc:0.6f}_epoch{epoch+1}'
            save_path = Path(rundir) / file_name
            torch.save(model.state_dict(), save_path)
                      
    
if __name__ == '__main__':
        
    np.random.seed(1)
    torch.manual_seed(1)
    
    train(rundir='models', path='data', epochs=300, learning_rate=0.001, use_gpu=False)
