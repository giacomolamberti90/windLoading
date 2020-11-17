import numpy as np
import pickle
import torch
import time

from datetime import datetime
from pathlib import Path

from evaluate import run_model
from loader import load_data_1D, load_data_2D, load_data_tiles, face
from model import RMSNet, RMSCNNNet

from torch.utils.data import DataLoader
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score, mean_squared_error, mean_squared_log_error
from sklearn.utils import resample

from utils import *
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE

def train_NN(rundir, path, epochs, learning_rate, use_gpu):
    
    data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_1D(path)
    #data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_2D(path, 'train')
    #data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_tiles(path)
    
    train_loader = DataLoader(data_train, batch_size=len(data_train), num_workers=1, shuffle=True)
    valid_loader = DataLoader(data_valid, batch_size=len(data_valid), num_workers=1, shuffle=True)
    
    model = RMSNet()
    
    if use_gpu:
        model = model.cuda()

    optimizer = torch.optim.Adam(model.parameters(), learning_rate, weight_decay=0.0)
    #scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=.3, threshold=1e-5)

    best_val_mse = float('inf')

    start_time = datetime.now()
    
    for epoch in range(epochs):
                
        change = datetime.now() - start_time
        print('starting epoch {}. time passed: {}\n'.format(epoch+1, str(change)))
        
        train_loss, train_mse, _, _ = run_model(model, train_loader, train=True, optimizer=optimizer)
        print(f'train loss: {train_loss:0.6f}')
        print(f'train MSE: {train_mse:0.6f}\n')         
        
        val_loss, val_mse, _, _ = run_model(model, valid_loader)
        print(f'valid loss: {val_loss:0.6f}')
        print(f'valid MSE: {val_mse:0.6f}\n')
        
        #scheduler.step(val_loss)
        
        if val_mse < best_val_mse:
            best_val_mse = val_mse

            file_name = f'val{val_mse:0.6f}_train{train_mse:0.6f}_epoch{epoch+1}'
            save_path = Path(rundir) / file_name
            torch.save(model.state_dict(), save_path)
                        

def train(path, use_gpu, method):
    
    data_train, data_valid, data_00deg, data_20deg, data_40deg, data_60deg, data_80deg = load_data_1D(path)
    #data_train, data_valid, data_00deg, data_20deg, data_45deg = load_data_tiles(path)
        
    if method == "linear":
        model = LinearRegression()
        
    if method == "quadratic":            
        polynomial_features= PolynomialFeatures(degree=2)    
        data_train.features = polynomial_features.fit_transform(data_train.features)
        data_valid.features = polynomial_features.fit_transform(data_valid.features)
        model = LinearRegression()#Ridge(alpha=0.01)
        
    if method == "ridge":
        model = Ridge(alpha=0.01)
        
    if method == "lasso":
        model = Lasso(alpha=0.01)
        
    if method == "randomForest":
        model = RandomForestRegressor(n_estimators=1000, max_features=2, max_depth=None, random_state=10)
    
    if method == "boosting":
        model = GradientBoostingRegressor(loss='quantile', n_estimators=100, max_depth=None, random_state=1)
    
    if method == "neuralNet":
        model = MLPRegressor(hidden_layer_sizes=(10,) * 5, activation='relu', solver='adam', 
                         alpha=0.01, learning_rate='adaptive', learning_rate_init=0.001)        
        #model.out_activation_ = 'relu'
        
    print('Total number of features: %d' % data_train.features.shape[1])
    
    # boot-strap
    Nbootstrap = 1
    coeff = []
    allImportance = []
    
    allPreds_train = []
    allPreds_valid = []
    
    allPreds_00deg = []
    allPreds_20deg = []
    allPreds_40deg = []
    allPreds_60deg = []
    allPreds_80deg = []
    
    for i in range(Nbootstrap):
        print(f'Bootstrap sample: {i:6d}')
        t = time.time()
        N = len(data_train) // len(data_00deg)
        train = resample(np.arange(0, len(data_train)), n_samples=len(data_train) // N)
        model.fit(data_train.features[train], data_train.labels[train])
        #model.fit(data_train.features, data_train.labels)
        allPreds_train.append(model.predict(data_train.features))
        allPreds_valid.append(model.predict(data_valid.features))

        #allPreds_00deg.append(model.predict(data_00deg.features))
        #allPreds_20deg.append(model.predict(data_20deg.features))
        #allPreds_40deg.append(model.predict(data_40deg.features))
        #allPreds_60deg.append(model.predict(data_60deg.features))
        allPreds_80deg.append(model.predict(data_80deg.features))
        
        elapsed = time.time() - t
        print(f"Time: {elapsed:6.4f}\n")
        #allImportance.append(garson(model.coefs_[0], model.coefs_[-1]))
    
    # save data
    #np.savetxt("80deg/prediction_80deg_train90.txt", allPreds_80deg)
    
    """
    allPreds_00deg = np.loadtxt("00deg/prediction_00deg_train10+90_cf.txt")
    allPreds_20deg = np.loadtxt("20deg/prediction_20deg_train10+30_cf.txt")
    allPreds_40deg = np.loadtxt("40deg/prediction_40deg_train30+50_cf.txt")
    allPreds_60deg = np.loadtxt("60deg/prediction_60deg_train50_cf.txt")
    """
    allPreds_80deg_lowestDKL = np.loadtxt("80deg/prediction_80deg_train70+90.txt")
    
    # train -------------------------------------------------------------------
    preds_mean = np.mean(np.asarray(allPreds_train), axis=0)
    
    labels = data_train.labels
    mse = np.sqrt(np.mean((labels-preds_mean)**2))
    
    print(f'RMSE train: {mse:0.6f}')

    preds = data_train.empiricalModels[:,0]
    mse = np.sqrt(np.mean((labels-preds)**2))
    print(f'RMSE train Paterson: {mse:0.6f}\n')

    # valid -------------------------------------------------------------------
    #postProcessing(np.asarray(allPreds_valid), data_valid, 45, method)
    
    # test performance --------------------------------------------------------    
    #postProcessing(np.asarray(allPreds_00deg), np.asarray(allPreds_00deg), data_00deg, 0, method)
    #postProcessing(np.asarray(allPreds_20deg), np.asarray(allPreds_20deg), data_20deg, 20, method)
    #postProcessing(np.asarray(allPreds_40deg), np.asarray(allPreds_40deg), data_40deg, 40, method)
    #postProcessing(np.asarray(allPreds_60deg), np.asarray(allPreds_60deg), data_60deg, 60, method)
    postProcessing(np.asarray(allPreds_80deg_lowestDKL), np.asarray(allPreds_80deg), data_80deg, 80, method)
    
if __name__ == '__main__':
        
    np.random.seed(1)
    torch.manual_seed(1)
    
    #train_NN(rundir='models', path='data', epochs=100, learning_rate=0.1, use_gpu=False)
    train(path='data', use_gpu=False, method="neuralNet")
