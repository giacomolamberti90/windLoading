import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker

from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels

def plot_class(split, angle, face, preds, labels, filename):
    
    panels = np.loadtxt('data/windTunnel/designPressure_' + split)    
    indices = np.where((panels[:,4] == angle) & (panels[:,5] == face))[0]
    
    accuracy = np.mean(preds[indices] == labels[indices])
    
    # make contour plot
    plt.rcParams.update({'font.size': 15})   
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True)
    
    ax = axs[0]
    # make contour plot    
    for i in indices:
                
        x = [panels[i,0], panels[i,1], panels[i,1], panels[i,0]]
        y = [panels[i,2], panels[i,2], panels[i,3], panels[i,3]]
        
        if preds[i] == 0:
            ax.fill(x,y,'g',edgecolor='k')
        elif preds[i]  == 1:
            ax.fill(x,y,'y',edgecolor='k')
        elif preds[i]  == 2:
            ax.fill(x,y,'r',edgecolor='k')
        elif preds[i]  == 3:
            ax.fill(x,y,'k',edgecolor='k')

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 2])
        
        ax.set_xlabel(r"$x[m]$"); ax.set_ylabel(r"$y[m]$")
        ax.set_title(r"$CNN$")
    
    ax = axs[1]
    # make contour plot    
    for i in indices:
        
        x = [panels[i,0], panels[i,1], panels[i,1], panels[i,0]]
        y = [panels[i,2], panels[i,2], panels[i,3], panels[i,3]]
        
        if labels[i] == 0:
            ax.fill(x,y,'g',edgecolor='k')
        elif labels[i] == 1:
            ax.fill(x,y,'y',edgecolor='k')
        elif labels[i] == 2:
            ax.fill(x,y,'r',edgecolor='k')
        elif labels[i] == 3:
            ax.fill(x,y,'k',edgecolor='k')

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 2])
        
        ax.set_xlabel(r"$x[m]$")
        ax.set_title(r"$PoliMi$")

    ax.text(0.05, 0.5, "accuracy:%6.2f" % accuracy)
    
    #ax.text(0.05, 0.3, "Test  acc:%6.2f" % accuracy[1])
    plt.savefig(filename)
    
    return accuracy
        
def plot_panel(images, feature, tile):
        
    # make contour plot
    plt.rcParams.update({'font.size': 18})   
    plt.figure(figsize=(3,5))

    cmin = np.min(np.asarray(images)[:,:,:,feature])
    cmax = np.max(np.asarray(images)[:,:,:,feature])
    
    print(cmin)
    print(cmax)
    
    for image in images:
        if tile == 'D':
            x = image[:,:,2]
        else:
            x = image[:,:,0]
        
        y = image[:,:,1]
        z = image[:,:,feature]
        
        plt.contourf(x, y, z, cmap='hot', edgecolor='k')
        plt.clim([cmin, cmax])
        
    #plt.xlim([0, 1])
    #plt.ylim([0, 2])
    
def plot_confusion_matrix(y_true, y_pred, classes, filename,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    plt.savefig(filename)    
    
    return ax    