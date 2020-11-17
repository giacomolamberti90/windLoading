import numpy as np
import pandas as pd
#import active_subspaces as ac

from misc import process_inputs_outputs, process_inputs
from loader import load_data

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from utils import *

def compute(X, f, weights=None):
        # Check inputs
        if X is not None:
            X, M, m = process_inputs(X)
            
        if weights is None:
            # default weights is for Monte Carlo
            weights = np.ones((M, 1)) / M
                             
        if X is None or f is None:
            raise Exception('X or f is None')
        e, W = ols_subspace(X, f, weights)

        return e, W                     
                             
def sorted_eigh(C):
    """Compute eigenpairs and sort.
    
    Parameters
    ----------
    C : ndarray
        matrix whose eigenpairs you want
        
    Returns
    -------
    e : ndarray
        vector of sorted eigenvalues
    W : ndarray
        orthogonal matrix of corresponding eigenvectors
    
    Notes
    -----
    Eigenvectors are unique up to a sign. We make the choice to normalize the
    eigenvectors so that the first component of each eigenvector is positive.
    This normalization is very helpful for the bootstrapping. 
    """
    e, W = np.linalg.eigh(C)
    e = abs(e)
    ind = np.argsort(e)
    e = e[ind[::-1]]
    W = W[:,ind[::-1]]
    s = np.sign(W[0,:])
    s[s==0] = 1
    W = W*s
    
    return e.reshape((e.size,1)), W

def ols_subspace(X, f, weights):
    """Estimate one-dimensional subspace with global linear model.
    
    Parameters
    ----------
    X : ndarray
        M-by-m matrix of input samples, oriented as rows
    f : ndarray
        M-by-1 vector of output samples corresponding to the rows of `X`
    weights : ndarray
        M-by-1 weight vector, corresponds to numerical quadrature rule used to
        estimate matrix whose eigenspaces define the active subspace
        
    Returns
    -------
    e : ndarray
        m-by-1 vector of eigenvalues
    W : ndarray
        m-by-m orthogonal matrix of eigenvectors
        
    Notes
    -----
    Although the method returns a full set of eigenpairs (to be consistent with
    the other subspace functions), only the first eigenvalue will be nonzero,
    and only the first eigenvector will have any relationship to the input 
    parameters. The remaining m-1 eigenvectors are only orthogonal to the first.
    """
    X, f, M, m = process_inputs_outputs(X, f)
    
    # solve weighted least squares
    A = np.hstack((np.ones((M, 1)), X)) * np.sqrt(weights)
    b = f * np.sqrt(weights)
    u = np.linalg.lstsq(A, b)[0]
    w = u[1:].reshape((m, 1))
    
    # compute rank-1 C
    C = np.dot(w, w.transpose())
    
    return sorted_eigh(C)


if __name__ == "__main__":
    
    dataFrame_00deg, features_00deg, labels_00deg = load_data(0)
    dataFrame_10deg, features_10deg, labels_10deg = load_data(10)
    dataFrame_30deg, features_30deg, labels_30deg = load_data(30)
    dataFrame_40deg, features_40deg, labels_40deg = load_data(40)
    dataFrame_90deg, features_90deg, labels_90deg = load_data(90)
    
    # PCA
    X_00deg = PCA(n_components=2).fit_transform(features_00deg)
    X_10deg = PCA(n_components=2).fit_transform(features_10deg)
    X_30deg = PCA(n_components=2).fit_transform(features_30deg)
    X_40deg = PCA(n_components=2).fit_transform(features_40deg)    
    X_90deg = PCA(n_components=2).fit_transform(features_90deg)
    
    # TSNE
    #X_00deg = TSNE(n_components=2).fit_transform(features_00deg)
    #X_10deg = TSNE(n_components=2).fit_transform(features_10deg)
    #X_90deg = TSNE(n_components=2).fit_transform(features_90deg)
    
    plt.rcParams.update({'font.size': 18})
    plt.figure(figsize=(5,5))
    plt.plot(X_00deg[:,0], X_00deg[:,1], '.r', markersize=12)
    plt.plot(X_10deg[:,0], X_10deg[:,1], '.g', markersize=12)
    plt.plot(X_90deg[:,0], X_90deg[:,1], '.b', markersize=12)
    plt.legend([r"$0^\circ$", r"$10^\circ$", r"$90^\circ$"])
    plt.xlabel(r"$X_1$")
    plt.ylabel(r"$X_2$")
    #plt.savefig("TSNE_2comp_all_30+40.png")
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(X_00deg[0::10,0], X_00deg[0::10,1], labels_00deg[0::10], c='r', marker='.')
    ax.scatter(X_10deg[0::10,0], X_10deg[0::10,1], labels_10deg[0::10], c='g', marker='.')
    ax.scatter(X_90deg[0::10,0], X_90deg[0::10,1], labels_90deg[0::10], c='b', marker='.')

    ax.set_xlabel(r"$X_1$")
    ax.set_ylabel(r"$X_2$")
    ax.set_zlabel(r"$C_p'$")
    
    plt.legend([r"$0^\circ$", r"$10^\circ$", r"$90^\circ$"])
    plt.show()
    
    '''
    N, m = features.shape
    NN = N
    
    ind = np.random.randint(low=0, high=N, size=(NN,))
    
    features = features[ind,:]
    labels = labels[ind]
    labels = labels.reshape(NN, 1)
    
    # global linear models with expectations and standard dev
    eigenvals, eigenvecs = compute(features, labels)
    w = eigenvecs[:,0].reshape((m, 1))
    
    plt.figure(figsize=(5,5))
    plt.plot(range(1, m+1), w, 'ko-', markersize=12)
    plt.xlabel('Weights')
    plt.ylabel('N')
    plt.grid(True)
    plt.axis([1, m, -1, 1])
    #plt.savefig('weights_' + str(ang) + 'deg_face' + str(face) + '_allFeatures.png')
    
    y = np.dot(features, w)
    plt.figure(figsize=(5,5))
    plt.plot(y, labels, 'bo', markersize=12)
    plt.xlabel('Active variable')
    plt.ylabel('Labels')
    plt.grid(True)
    #plt.savefig('activeVariable_' + str(ang) + 'deg_face' + str(face)  + '_allFeatures.png')
    
    sv = np.linalg.svd(features, full_matrices=False, compute_uv=False)
    print(sv)
    '''