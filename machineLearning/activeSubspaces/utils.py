import numpy as np
import scipy

import matplotlib.pyplot as plt
from matplotlib import ticker

import scipy.stats as stats

def plot_contour(coords, variable, clim, filename):
        
    # make contour plot
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
 
def plot_profiles(data, preds, se, plane, angle):
    
    plt.rcParams.update({'font.size': 18})       
    
    coords = data.coords
    labels = data.labels
    faces = data.faces
    paterson_holmes = data.empiricalModels[:,0]
    #selvam = data.empiricalModels[:,1]
    #paterson = data.empiricalModels[:,2]
    #richards = data.empiricalModels[:,3]
    
    X = np.zeros_like(faces)
    Y = np.zeros_like(faces)
    
    if plane == 'Ynorm':        
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

        Y_values = [0.503, 0.997, 1.484]

    if plane == 'Znorm':
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
        
        idx = np.round(Y,3) == Y_values[i]
        
        X_idx = X[idx]
        labels_idx = labels[idx]
        preds_idx = preds[idx]
        se_idx = se[idx]
        paterson_holmes_idx = paterson_holmes[idx]
        
        idx = np.argsort(X_idx)
        
        plt.figure()
        plt.plot(X_idx[idx], labels_idx[idx], '.r')
        plt.plot(X_idx[idx], preds_idx[idx], 'k')
        
        plt.fill_between(X_idx[idx], preds_idx[idx] - 2*se_idx[idx], preds_idx[idx] + 2*se_idx[idx], 
                          color='k', alpha=.2)
        
        plt.plot(X_idx[idx], paterson_holmes_idx[idx], '.b')     
        plt.xlabel(r"$x[m]$"); plt.ylabel(r"$C_p'$")
        #if i == 0 and plane != 'Znorm':
        #    plt.legend(['LES', 'Linear-regression', 'Paterson-Holmes'])            
        plt.tight_layout()

        #plt.title(filename)
        if plane == 'Znorm':
            plt.title(r"$Z=%6.1f$" % Y_values[i])
            plt.ylim([0, 0.4])
            plt.savefig('profile_3_neuralNet_' + str(angle) + 'deg.png')
        else:
            plt.title(r"$Y=%6.1f$" % Y_values[i])
            plt.ylim([0, 0.3])
            plt.savefig('profile_' + str(i) + '_neuralNet_' + str(angle) + 'deg.png')

def kl_divergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))
    
def plotDensity(x, y, angle):

    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 25})
             
    xmin = -5
    xmax = 2
    ymin = -1.5
    ymax = 1

    # Create meshgrid
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = stats.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    
    plt.figure()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    c = plt.contourf(xx, yy, f, cmap='coolwarm')
    #cbar = plt.colorbar(c)
    #cset = plt.contour(xx, yy, f, colors='k')
    #plt.clabel(cset, inline=1, fontsize=18)
    plt.xlabel('$X_1$')
    plt.ylabel('$X_2$')
    plt.tight_layout()
    plt.savefig('features_pdf_' + str(angle) + 'deg.png')

    return f

    