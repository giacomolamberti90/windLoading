import numpy as np
import scipy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker

from scipy.interpolate import griddata

def plot_contour(variable, skip, clim, filename, title):
    
    coords = np.loadtxt('coords.txt')[::skip,:]
    size = 100
    
    idx = variable != 0
    
    variable = variable[idx]
    coords = coords[idx]
    
    # make contour plot
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 18})   
    plt.figure(figsize=(7,5))
    
    # left
    index = np.nonzero(coords[:,2] <= 1e-3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]
    plt.scatter(-x-0.2, y, s=size, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(-0.8, -0.2, r"$left$")    
    plt.xlabel(r"$x[m]$")
        
    # front
    index = np.nonzero(coords[:,0] <= 1e-4)
    x = coords[index,2]
    y = coords[index,1] 
    z = variable[index]
    plt.scatter(x, y, s=size, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(-0.1, -0.2, r"$front$")
    plt.ylabel(r"$y[m]$")
    
    # top
    index = np.nonzero(coords[:,1] == 2)
    x = coords[index,2]
    y = coords[index,0]        
    z = variable[index]
    plt.scatter(x, y+2.2, s=size, c=z, marker='.', cmap='hot')            
    plt.clim(clim)
    plt.text(0, 3.3, r"$top$")
    
    # right
    index = np.nonzero(coords[:,2] >= 0.3)
    x = coords[index,0]
    y = coords[index,1]
    z = variable[index]    
    plt.scatter(x+0.5, y, s=size, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(0.8, -0.2, r"$right$")
    plt.xlabel(r"$x[m]$")
    
    # back
    index = np.nonzero(coords[:,0] >= 1-1e-4)
    x = coords[index,2]
    y = coords[index,1]
    z = variable[index]
    plt.scatter(-x+2, y, s=size, c=z, marker='.', cmap='hot')
    plt.clim(clim)
    plt.text(1.6, -0.2, r"$back$")
    
    cb = plt.colorbar(); cb.locator = ticker.MaxNLocator(nbins=5); cb.update_ticks()
    plt.xticks([])
    plt.yticks([])
    plt.clim(clim)
    #plt.title(title)
    #plt.axis('equal')
    plt.tight_layout()
    ax = plt.gca()
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.show()
    #plt.savefig(filename)
    
def plot_contourf(variable, skip, clim, filename, title):
    
    coords = np.loadtxt('coords.txt')[::skip,:]
    levels = 100
    
    idx = variable != 0
    
    variable = variable[idx]
    coords = coords[idx]
    
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
    
def plot_profiles(preds, low, up, skip, plane, color, ylim):
    
    plt.rcParams.update({'font.size': 18})
    
    edges = np.loadtxt('edges.txt')
    coords = np.loadtxt('coords.txt')[edges == 0]
    faces = np.loadtxt('faces.txt')[edges == 0]

    preds = griddata(np.loadtxt('coords.txt')[::skip,:], preds, coords, method='linear')
    
    X = np.zeros_like(faces)
    Y = np.zeros_like(faces)
    
    if plane == '0':        
        for i, face in enumerate(faces):
            
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
        for i, face in enumerate(faces):
            
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
        for i, face in enumerate(faces):
            
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
        for i, face in enumerate(faces):
                   
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
                            
            preds_idx = preds[idx]            
            idx = np.argsort(X_idx)
            
            #plt.figure()
            if color[0] == '-':
                lw = 2
            else:
                lw = 1
                
            plt.plot(X_idx[idx], preds_idx[idx], color, linewidth=lw)
            
            plt.ylim(ylim)
            
            plt.xlabel(r"$x[m]$"); plt.ylabel(r"$C_P$")
            plt.tight_layout()
            