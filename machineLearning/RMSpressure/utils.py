import numpy as np
import scipy


import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker

def plot_contour(coords, variable, clim, filename):
        
    # make contour plot
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 18})   
    plt.figure(figsize=(7,5))
    
    # left
    index = np.nonzero(coords[:,2] <= 1e-3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]
    plt.scatter(-x-0.2, y, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(-0.8, -0.2, r"$left$")    
    plt.xlabel(r"$x[m]$")
        
    # front
    index = np.nonzero(coords[:,0] <= 1e-4)
    x = coords[index,2]
    y = coords[index,1] 
    z = variable[index]    
    plt.scatter(x, y, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(-0.1, -0.2, r"$front$")
    plt.ylabel(r"$y[m]$")
    
    # top
    index = np.nonzero(coords[:,1] == 2)
    x = coords[index,2]
    y = coords[index,0]        
    z = variable[index]
    plt.scatter(x, y+2.2, c=z, marker='.', cmap='hot')            
    plt.clim(clim)
    plt.text(0, 3.3, r"$top$")
    
    # right
    index = np.nonzero(coords[:,2] >= 0.3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]    
    plt.scatter(x+0.5, y, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(0.8, -0.2, r"$right$")
    plt.xlabel(r"$x[m]$")
    
    # back
    index = np.nonzero(coords[:,0] >= 1-1e-4)
    x = coords[index,2]
    y = coords[index,1]
    z = variable[index]
    plt.scatter(-x+2, y, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(1.6, -0.2, r"$back$")
    
    cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    plt.xticks([])
    plt.yticks([])
    plt.clim(clim)
    #plt.title(filename)
    #plt.axis('equal')
    plt.tight_layout()
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    #plt.show()
    plt.savefig(filename)
    

def plot_contourf(coords, variable, clim, filename):
    
    levels = 100
    
    #idx = variable > 0
    
    #variable = variable[idx]
    #coords = coords[idx]
    
    # make contour plot
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 18})   
    plt.figure(figsize=(7,5))
    
    # left
    index = np.nonzero(coords[:,2] <= 1e-3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(x, y, z)
    zi = rbf(xi,yi)
    plt.contourf(-xi-0.2, yi, zi, levels, cmap='hot')
    plt.clim(clim)
    plt.text(-0.8, -0.2, r"$left$")    
        
    # front
    index = np.nonzero(coords[:,0] <= 1e-4)
    x = coords[index,2]
    y = coords[index,1] 
    z = variable[index]
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(x, y, z)
    zi = rbf(xi,yi)
    plt.contourf(xi, yi, zi, levels, cmap='hot')    
    plt.clim(clim)
    plt.text(-0.1, -0.2, r"$front$")
    plt.ylabel(r"$y[m]$")
    
    # top
    index = np.nonzero(coords[:,1] == 2)
    x = coords[index,2]
    y = coords[index,0]    
    z = variable[index]
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(x, y, z)
    zi = rbf(xi,yi)
    plt.contourf(xi, yi+2.2, zi, levels, cmap='hot')    
    plt.clim(clim)
    plt.text(0, 3.3, r"$top$")
    
    # right
    index = np.nonzero(coords[:,2] >= 0.3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(x, y, z)
    zi = rbf(xi,yi)
    plt.contourf(xi+0.5, yi, zi, levels, cmap='hot')
    plt.clim(clim)
    plt.text(0.8, -0.2, r"$right$")
    plt.text(0.2, -0.4, r"$x[m]$")
    #plt.xlabel(r"$x[m]$")
    
    # back
    index = np.nonzero(coords[:,0] >= 1-1e-4)
    x = coords[index,2]
    y = coords[index,1]
    z = variable[index]
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)
    rbf = scipy.interpolate.Rbf(x, y, z)
    zi = rbf(xi,yi)
    plt.contourf(-xi+2, yi, zi, levels, cmap='hot')    
    plt.clim(clim)
    plt.text(1.6, -0.2, r"$back$")
    
    plt.xticks([])
    plt.yticks([])
    plt.clim(clim)    
    #cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    #plt.title(title)
    #plt.axis('equal')
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.tight_layout()

    ax, _ = mpl.colorbar.make_axes(plt.gca())
    cbar = mpl.colorbar.ColorbarBase(ax, cmap='hot', norm=mpl.colors.Normalize(vmin=clim[0], vmax=clim[1]))
    cbar.set_clim(clim)

    #plt.show()
    plt.savefig(filename, bbox_inches = "tight")    
    
    
def plot_contour_face(coords, variable, clim, face, filename):

    plt.rcParams.update({'font.size': 18})       
    
    if face == 1:
        # front
        plt.figure(figsize=(3,5))
        index = np.nonzero(coords[:,0] <= 1e-4)
        x = coords[index,2]
        y = coords[index,1] 
        z = variable[index]    
        plt.scatter(x, y, c=z, marker='.', cmap='hot')
        plt.clim(clim)
        plt.ylabel(r"$y[m]$")
    
    if face == 2:
        # top
        plt.figure(figsize=(3,4))
        index = np.nonzero(coords[:,1] == 2)
        x = coords[index,2]
        y = coords[index,0]        
        z = variable[index]
        plt.scatter(x, y, c=z, marker='.', cmap='hot')            
        plt.clim(clim)
        
    if face == 3:
        # left
        plt.figure(figsize=(4,5))
        index = np.nonzero(coords[:,2] <= 1e-3)
        x = coords[index,0]
        y = coords[index,1]
        z = variable[index]
        plt.scatter(-x, y, c=z, marker='.', cmap='hot')
        plt.clim(clim)
        plt.xlabel(r"$x[m]$")
    
    if face == 4:
        # back
        plt.figure(figsize=(3,5))
        index = np.nonzero(coords[:,0] >= 1-1e-4)
        x = coords[index,2]
        y = coords[index,1]
        z = variable[index]
        plt.scatter(-x, y, c=z, marker='.', cmap='hot')
        plt.clim(clim)
        
    if face == 5:
        # right
        plt.figure(figsize=(4,5))
        index = np.nonzero(coords[:,2] >= 0.3)
        x = coords[index,0]
        y = coords[index,1]
        z = variable[index]    
        plt.scatter(x, y, c=z, marker='.', cmap='hot')
        plt.clim(clim)
        plt.xlabel(r"$x[m]$")
        
    cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    plt.xticks([])
    plt.yticks([])
    plt.clim(clim)
    #plt.title(filename)
    #plt.axis('equal')
    plt.tight_layout()
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.show()
    plt.savefig(filename)
    
def plot_profiles(data, preds, low, up, plane, color, alpha):
    
    plt.rcParams.update({'font.size': 18})
    
    coords = data.coords
    labels = data.labels
    faces = data.faces
    paterson_holmes = data.empiricalModels[:,0]

    X = np.zeros_like(faces)
    Y = np.zeros_like(faces)
    
    if plane == '0':        
        for i,face in enumerate(faces):
            
            if face == 1:
                X[i] = coords[i,2]
                Y[i] = coords[i,1]
            if face == 5:
                X[i] = coords[i,0] + 0.3
                Y[i] = coords[i,1]
            if face == 4:
                X[i] = coords[i,2] + 1.3
                Y[i] = coords[i,1]
            if face == 3:
                X[i] = coords[i,0] + 1.6
                Y[i] = coords[i,1]

        Y_values = [0.503]
    
    if plane == '1':        
        for i,face in enumerate(faces):
            
            if face == 1:
                X[i] = coords[i,2]
                Y[i] = coords[i,1]
            if face == 5:
                X[i] = coords[i,0] + 0.3
                Y[i] = coords[i,1]
            if face == 4:
                X[i] = coords[i,2] + 1.3
                Y[i] = coords[i,1]
            if face == 3:
                X[i] = coords[i,0] + 1.6
                Y[i] = coords[i,1]

        Y_values = [0.997]
    
    if plane == '2':        
        for i,face in enumerate(faces):
            
            if face == 1:
                X[i] = coords[i,2]
                Y[i] = coords[i,1]
            if face == 5:
                X[i] = coords[i,0] + 0.3
                Y[i] = coords[i,1]
            if face == 4:
                X[i] = coords[i,2] + 1.3
                Y[i] = coords[i,1]
            if face == 3:
                X[i] = coords[i,0] + 1.6
                Y[i] = coords[i,1]

        Y_values = [1.484]
        
    if plane == '3':
        for i,face in enumerate(faces):
                   
            if face == 1:
                X[i] = coords[i,1]
                Y[i] = coords[i,2]
            if face == 2:
                X[i] = coords[i,0] + 2
                Y[i] = coords[i,2]
            if face == 4:
                X[i] = coords[i,1] + 3
                Y[i] = coords[i,2]
                
        Y_values = [0.156]
             
    for i in range(len(Y_values)):
        
        for face in range(6):
            idx = (np.round(Y,3) == Y_values[i]) & (faces == face)
            
            if face == 3: # left face
                X_idx = 4.2 - X[idx]
            else:
                X_idx = X[idx]
            
            labels_idx = labels[idx]
            preds_idx = preds[idx]
            low_idx = low[idx]
            up_idx = up[idx]
            paterson_holmes_idx = paterson_holmes[idx]
            
            idx = np.argsort(X_idx)
            
            plt.plot(X_idx[idx], labels_idx[idx], '.r')
            plt.plot(X_idx[idx], preds_idx[idx], color)
            
            plt.fill_between(X_idx[idx], low_idx[idx], up_idx[idx], color=color, alpha=alpha)
            
            plt.plot(X_idx[idx], paterson_holmes_idx[idx], '.b')
            
            plt.xlabel(r"$x[m]$"); plt.ylabel(r"$C_p'$")
            #if i == 0 and plane != 'Znorm':
            #    plt.legend(['LES', 'Linear-regression', 'Paterson-Holmes'])            
            plt.tight_layout()
            plt.ylim([-0.01, 0.4])
    
def plot(labels, preds, filename):

    mse = np.sqrt(np.mean((labels-preds)**2))
    
    plt.rcParams.update({'font.size': 18})

    plt.figure(figsize=(5,5))
    plt.plot(labels, preds, '.r')
    plt.plot(np.linspace(0,0.5), np.linspace(0,0.5), '-k')
    plt.plot(np.linspace(0,0.5), 1.1*np.linspace(0,0.5), '--k')
    plt.plot(np.linspace(0,0.5), 0.9*np.linspace(0,0.5), '--k')
    plt.xlabel(r"$C_p'$ - LES"); plt.ylabel(r"$C_p'$ - pred")
    plt.axis([0, 0.3, 0, 0.3])
    plt.text(0.025, 0.275, "RMSE = %6.4f" % mse)
    plt.tight_layout()
    plt.savefig(filename)

def garson(A, B):
    """
    Computes Garson's algorithm
    A = matrix of weights of input-hidden layer (rows=input & cols=hidden)
    B = vector of weights of hidden-output layer
    """
    #B = np.diag(B)

    # connection weight through the different hidden node
    cw = np.dot(A, B)

    # weight through node (axis=0 is column; sum per input feature)
    cw_h = abs(cw).sum(axis=0)

    # relative contribution of input neuron to outgoing signal of each hidden neuron
    # sum to find relative contribution of input neuron
    rc = np.divide(abs(cw), abs(cw_h))
    rc = rc.sum(axis=1)

    # normalize to 100% for relative importance
    ri = rc / rc.sum()
    
    return(ri)


def performance(pred_mean, pred_up, pred_low, y):
    
    accuracy = 0
    err = 0
    for i in range(len(y)):
    
        if y[i] > pred_up[i]:
            if y[i] - pred_up[i] > err:
                err = y[i] - pred_up[i]
            
        elif y[i] < pred_low[i]:
            if pred_low[i] - y[i] > err:
                err = pred_low[i] - y[i]
                
        else:
            accuracy += 1      
    
    return err, accuracy/len(y)


def postProcessing(allPreds, allPreds_best, data, angle, method):
    
    preds_mean = np.mean(allPreds, axis=0)    
    preds_se = np.std(allPreds, axis=0)

    preds_low = np.percentile(allPreds, 2.5, axis=0)
    preds_up = np.percentile(allPreds, 97.5, axis=0)
    
    preds_mean_best = np.mean(allPreds_best, axis=0)    
    preds_se_best = np.std(allPreds_best, axis=0)

    preds_low_best = np.percentile(allPreds_best, 2.5, axis=0)
    preds_up_best = np.percentile(allPreds_best, 97.5, axis=0)    
    
    '''
    plt.figure()
    variables = [r"$C_P$", r"$\frac{k}{U_{ref}^2}$", r"$\frac{U_0}{U_{ref}^2}$", 
                 r"$\frac{H \vert\vert\nabla P\vert\vert}{0.5\rho U_{ref}^2}$", r"$C_f$"]
    #variables = [r"$C_P$", r"$k$", r"$U_0$", r"$y^+$", r"$T$", r"$\nu_t$", r"$\vert\vert\nabla P\vert\vert$", r"$C_f$"]
    plt.bar(range(len(variables)), importance, color=['black', 'red', 'blue', 'green', 'yellow'])
    #plt.bar(range(len(variables)), importance)
    plt.xticks(range(len(variables)), variables)
    plt.tight_layout()
    plt.savefig('importance_' + str(angle) + '.png')
 
    # plot coefficients -------------------------------------------------------
    m = data_train.features.shape[1]
    plt.figure(figsize=(5,5))
    plt.plot(range(1, m+1), coeff_mean, 'ko-', markersize=12)
    plt.xlabel('Weights')
    plt.ylabel('N')
    plt.grid(True)
    plt.axis([1, m, -1, 1])
    '''
    # test performance --------------------------------------------------------
    labels = data.labels
    mse = np.sqrt(np.mean((labels-preds_mean)**2))
    print(f'{angle:2d} RMSE test: {mse:0.6f}')

    preds_emp = data.empiricalModels[:, 0]
    mse = np.sqrt(np.mean((labels-preds_emp)**2))
    print(f'{angle:2d} RMSE test Paterson: {mse:0.6f}')

    err, _ = performance(preds_emp, preds_emp, preds_emp, labels)
    
    print(f'{angle:2d} Maximum discrepancy Paterson: {err:0.6f}')
    
    err, acc = performance(preds_mean, preds_up, preds_low, labels)
    
    print(f'{angle:2d} Maximum discrepancy: {err:0.6f}')
    print(f'{angle:2d} Points inside CI: {acc:0.6f}\n')
    
    #plot(labels, preds_mean, 'fit_' + method + '_' + str(angle) + 'deg.png')
    
    #plot_contour(data.coords, labels, [0, 0.3], '')
    plot_contour(data.coords, preds_mean, [0, 0.3], 'cp_rms_' + method + '_lowestDKL_' + str(angle) + 'deg.png')
    plot_contour(data.coords, preds_mean_best, [0, 0.3], 'cp_rms_' + method + '_best_' + str(angle) + 'deg.png')
    #plot_contourf(data.coords, preds_mean, [0, 0.3], 'cp_rms_' + method + '_' + str(angle) + 'deg.png')
    #plot_contour(data.coords, data.empiricalModels[:,0], [0, 0.3], '')
    
    
    plt.figure()
    plot_profiles(data, preds_mean, preds_low, preds_up, '0', 'Cyan', 0.4)
    plot_profiles(data, preds_mean_best, preds_low_best, preds_up_best, '0', 'k', 0.2)
    plt.savefig('profile_0_' + method + '_' + str(angle) + 'deg.png')
    
    plt.figure()
    plot_profiles(data, preds_mean, preds_low, preds_up, '1', 'Cyan', 0.4)
    plot_profiles(data, preds_mean_best, preds_low_best, preds_up_best, '1', 'k', 0.2)
    plt.savefig('profile_1_' + method + '_' + str(angle) + 'deg.png')
    
    plt.figure()
    plot_profiles(data, preds_mean, preds_low, preds_up, '2', 'Cyan', 0.4)
    plot_profiles(data, preds_mean_best, preds_low_best, preds_up_best, '2', 'k', 0.2)
    plt.savefig('profile_2_' + method + '_' + str(angle) + 'deg.png')
    
    plt.figure()
    plot_profiles(data, preds_mean, preds_low, preds_up, '3', 'Cyan', 0.4)
    plot_profiles(data, preds_mean_best, preds_low_best, preds_up_best, '3', 'k', 0.2)
    plt.savefig('profile_3_' + method + '_' + str(angle) + 'deg.png')