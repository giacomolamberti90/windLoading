from loader import load_data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from mpl_toolkits.mplot3d import Axes3D
from utils import *

def data_visualization(dataFrame_all):
    
    """
    Function that receives dataFrame and visualize data    
    """

    variables = ['cp_rms', 'cp_mean', 'tke', 'U_inflow', 'gradp', 'cf']
    
    summary = round(dataFrame_all[variables].describe(), 2)
    
    # histograms
    dataFrame_all[variables].hist(bins=15, color='steelblue', edgecolor='black', linewidth=1.0,
                             xlabelsize=8, ylabelsize=8, grid=False)
    plt.tight_layout(rect=(0, 0, 1.2, 1.2))
    
    # correlation matrix
    corr_matrix = dataFrame_all[variables].corr()
    
    f, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(round(corr_matrix, 2), annot=True, cmap="coolwarm",fmt='.2f',
                 linewidths=.05)
    f.suptitle('Correlation Matrix', fontsize=14)
    
    # pair-wise scatter plots
    sns.set(font_scale=2) 
    pp = sns.pairplot(data=dataFrame_all, vars=variables,hue='angle', size=1.8, aspect=1.8,
                      plot_kws=dict(edgecolor="k", linewidth=0.5),
                      diag_kind="kde", diag_kws=dict(shade=True))
    fig = pp.fig 
    fig.subplots_adjust(top=0.93, wspace=0.3)
    pp.savefig('pair-plot.png')
    
    # Scatter Plot
    plt.scatter(dataFrame_all['yPlus'], dataFrame_all['nut'],
                alpha=0.4, edgecolors='w')
    plt.xlabel('yPlus')
    plt.ylabel('nut')
    
def dimensionReduction(dataFrame, ncomp, skip, method):
    
    """
    Function that receives:
        - X: features matrix
        - ncomp: number of components in the reduced space
        - skip: number of datapoints to skip
        - method
        - return: new dataframe        
    """
    
    variables = ['cp_mean', 'tke', 'U_inflow', 'gradp', 'cf']
    X = dataFrame[variables].values
    #X = (X - np.mean(X, axis=0))/np.std(X, axis=0)
    y = dataFrame["cp_rms"].values
    faces = dataFrame["face"]
    
    if method == "PCA":
        pca = PCA(n_components=ncomp)
        X_reduced = pca.fit_transform(X[0::skip,])
        print(pca.explained_variance_ratio_)
    
    if method == "TSNE":
        tsne = TSNE(n_components=ncomp, verbose=1, perplexity=30, n_iter=300)
        X_reduced = tsne.fit_transform(X[0::skip,])
    
    y_reduced = y[0::skip].reshape(len(y[0::skip]), 1)
    f_reduced = faces[0::skip].reshape(len(y[0::skip]), 1)
    ang = dataFrame["angle"].values[0::skip].reshape(len(y[0::skip]), 1)
    
    dataFrame_reduced = pd.DataFrame(np.concatenate((y_reduced, X_reduced, ang, f_reduced), axis=1), 
                                     columns=['y', 'X1', 'X2', 'angle', 'face'])
    
    return dataFrame_reduced
    
def clusters(features_all, features_10deg, features_20deg, 
             features_30deg, features_45deg):
    
    """
    Split training data into clusters
    """
    
    # K_means
    model = KMeans(n_clusters=2, random_state=0)
    model.fit(features_all)
    
    # clusters
    labels_10deg = model.predict(features_10deg)
    labels_20deg = model.predict(features_20deg)
    labels_30deg = model.predict(features_30deg)
    labels_40deg = model.predict(features_40deg)
    
    print("Points in cluster 0 at 10deg:%6.2f" % (sum(labels_10deg == 0)/len(labels_10deg)))
    print("Points in cluster 0 at 20deg:%6.2f" % (sum(labels_20deg == 0)/len(labels_20deg)))
    print("Points in cluster 0 at 30deg:%6.2f" % (sum(labels_30deg == 0)/len(labels_30deg)))
    print("Points in cluster 0 at 45deg:%6.2f" % (sum(labels_40deg == 0)/len(labels_40deg)))
    
    print("Points in cluster 1 at 10deg:%6.2f" % (sum(labels_10deg == 1)/len(labels_10deg)))
    print("Points in cluster 1 at 20deg:%6.2f" % (sum(labels_20deg == 1)/len(labels_20deg)))
    print("Points in cluster 1 at 30deg:%6.2f" % (sum(labels_30deg == 1)/len(labels_30deg)))
    print("Points in cluster 1 at 45deg:%6.2f" % (sum(labels_40deg == 1)/len(labels_40deg)))
    
if __name__ == '__main__':

    # load data ===============================================================    
    dataFrame_00deg, features_00deg, labels_00deg = load_data(0)
    dataFrame_10deg, features_10deg, labels_10deg = load_data(10)
    dataFrame_20deg, features_20deg, labels_20deg = load_data(20)
    dataFrame_30deg, features_30deg, labels_30deg = load_data(30)
    dataFrame_40deg, features_40deg, labels_40deg = load_data(40)
    dataFrame_50deg, features_50deg, labels_50deg = load_data(50)
    dataFrame_60deg, features_60deg, labels_60deg = load_data(60)
    dataFrame_70deg, features_70deg, labels_70deg = load_data(70)
    dataFrame_80deg, features_80deg, labels_80deg = load_data(80)
    dataFrame_90deg, features_90deg, labels_90deg = load_data(90)
    
    dataFrame_list = [dataFrame_00deg, dataFrame_10deg, dataFrame_20deg, dataFrame_30deg, dataFrame_40deg,
                      dataFrame_50deg, dataFrame_60deg, dataFrame_70deg, dataFrame_80deg, dataFrame_90deg]
    dataFrame_all = pd.concat(dataFrame_list, ignore_index=True)
                    
    # dimension reduction =====================================================
    facestr = ""
    method = "PCA"
    
    dataFrame_DR = dimensionReduction(dataFrame_all, 2, 1, method)
    '''
    dataFrame_DR_00 = dimensionReduction(dataFrame_00deg, 2, 1, method)
    dataFrame_DR_10 = dimensionReduction(dataFrame_10deg, 2, 1, method)
    dataFrame_DR_30 = dimensionReduction(dataFrame_30deg, 2, 1, method)
    dataFrame_DR_40 = dimensionReduction(dataFrame_40deg, 2, 1, method)
    dataFrame_DR_90 = dimensionReduction(dataFrame_90deg, 2, 1, method)

    dataFrame_all = pd.concat([dataFrame_DR_00, dataFrame_DR_10, dataFrame_DR_30, dataFrame_DR_40, dataFrame_DR_90], ignore_index=True)
    '''
    dataFrame_DR_00 = dataFrame_DR[(dataFrame_DR['angle'] == 0) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_10 = dataFrame_DR[(dataFrame_DR['angle'] == 10) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_20 = dataFrame_DR[(dataFrame_DR['angle'] == 20) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_30 = dataFrame_DR[(dataFrame_DR['angle'] == 30) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_40 = dataFrame_DR[(dataFrame_DR['angle'] == 40) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_50 = dataFrame_DR[(dataFrame_DR['angle'] == 50) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_60 = dataFrame_DR[(dataFrame_DR['angle'] == 60) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_70 = dataFrame_DR[(dataFrame_DR['angle'] == 70) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_80 = dataFrame_DR[(dataFrame_DR['angle'] == 80) & (dataFrame_DR['face'] > 3)]
    dataFrame_DR_90 = dataFrame_DR[(dataFrame_DR['angle'] == 90) & (dataFrame_DR['face'] > 3)]
    
    # figures =================================================================
    #sns.set(font_scale=1.5)
    #dataFrame = dataFrame_DR[dataFrame_DR['angle'].isin([0, 10, 90])]    
    #ax = sns.lmplot(x="X1", y="X2", hue="angle", data=dataFrame)
    #plt.xlabel('$X_1$'); plt.ylabel('$X_2$')    
    #ax.savefig('lmplot_PCA_1_features+output.png')
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dataFrame_DR_00['X1'], dataFrame_DR_00['X2'], dataFrame_DR_00['y'], marker='.')
    ax.scatter(dataFrame_DR_10['X1'], dataFrame_DR_10['X2'], dataFrame_DR_10['y'], marker='.')
    ax.scatter(dataFrame_DR_90['X1'], dataFrame_DR_90['X2'], dataFrame_DR_90['y'], marker='.')
    plt.legend(['$0^\circ$', '$10^\circ$', '$90^\circ$'])
    ax.set_xlabel('$X_1$'); ax.set_ylabel('$X_2$'); ax.set_zlabel('$C_p''$')
    plt.savefig(method + '_features_3dplot_1' + facestr + '.png')
    '''
    angle = 80
    
    coords = dataFrame_00deg[['X', 'Y', 'Z']].values
    
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'font.size': 18})
    '''
    plot_contour(coords, dataFrame_DR_00['X1'].values, [-10, 10], 'X1_00deg.png')
    plot_contour(coords, dataFrame_DR_10['X1'].values, [-10, 10], 'X1_10deg.png')
    plot_contour(coords, dataFrame_DR_20['X1'].values, [-10, 10], 'X1_20deg.png')
    plot_contour(coords, dataFrame_DR_30['X1'].values, [-10, 10], 'X1_30deg.png')
    plot_contour(coords, dataFrame_DR_40['X1'].values, [-10, 10], 'X1_40deg.png')
    plot_contour(coords, dataFrame_DR_50['X1'].values, [-10, 10], 'X1_50deg.png')
    plot_contour(coords, dataFrame_DR_60['X1'].values, [-10, 10], 'X1_60deg.png')
    plot_contour(coords, dataFrame_DR_70['X1'].values, [-10, 10], 'X1_70deg.png')
    plot_contour(coords, dataFrame_DR_80['X1'].values, [-10, 10], 'X1_80deg.png')
    plot_contour(coords, dataFrame_DR_90['X1'].values, [-10, 10], 'X1_90deg.png')
    '''
    
    cells = np.zeros((len(dataFrame_DR_00), 1))
    cells[(dataFrame_DR_80['X1'] < 10) & (dataFrame_DR_80['X2'] < 0)] = 1
    plot_contour(coords, cells, [0, 1], 'highlight_80deg.png')
    
    plt.figure()
    plt.plot(dataFrame_DR_80['X1'], dataFrame_DR_80['X2'], '.')
    plt.plot(dataFrame_DR_70['X1'], dataFrame_DR_70['X2'], '.')
    plt.plot(dataFrame_DR_90['X1'], dataFrame_DR_90['X2'], '.')
    plt.legend(['$80^\circ$', '$70^\circ$', '$90^\circ$'])
    plt.xlabel("$X_1$"); plt.ylabel("$X_2$")
    plt.show()
    #plt.savefig(method + '_features_scatter_' + str(angle) + 'deg.png')
    
    plt.figure()
    dataFrame = dataFrame_DR[dataFrame_DR["angle"].isin([0, 30])]
    ax = sns.lmplot(x="X1", y="X2", hue="angle", data=dataFrame)
    plt.xlabel("$X_1$"); plt.ylabel("$X_2$")
    #plt.show()
    plt.savefig(method + '_features_lmplot_0+30deg.png')
    
    sns.set(font_scale=1.75)
    plt.figure()
    ax = sns.kdeplot(dataFrame_DR_00['X1'], dataFrame_DR_00['X2'], cmap="Reds", shade=True, shade_lowest=False)
    #ax = sns.kdeplot(dataFrame_DR_90['X1'], dataFrame_DR_90['y'], cmap="Reds", shade=True, shade_lowest=False) 
    #ax = sns.kdeplot(dataFrame_DR_70['X1'], dataFrame_DR_70['X2'], cmap="Blues", shade=True, shade_lowest=False) 
    #ax = sns.kdeplot(dataFrame_DR_80['X1'], dataFrame_DR_80['X2'], cmap="YlOrBr", shade=True, shade_lowest=False)
    plt.axis([-3, 3, -3, 3])
    #plt.legend(['$0^\circ$', '$10^\circ$', '$90^\circ$'])
    plt.xlabel("$X_1$"); plt.ylabel("$X_2$")
    plt.show()
    #plt.savefig(method + '_features_kde_' + str(angle) + 'deg.png')
    
    plt.figure()
    plt.hist(dataFrame_DR_00['X1'], bins='auto', alpha=0.8)
    plt.hist(dataFrame_DR_10['X1'], bins='auto', alpha=0.8)
    plt.hist(dataFrame_DR_90['X1'], bins='auto', alpha=0.8)
    plt.yscale('log')
    plt.legend(['$0^\circ$', '$10^\circ$', '$90^\circ$'])
    plt.xlabel("$X_1$")
    #plt.show()
    plt.savefig(method + '_features_log-hist_' + str(angle) + 'deg.png')
    
    plt.figure()
    ax = sns.kdeplot(dataFrame_DR_80['X1'], shade=True, color="blue")
    ax = sns.kdeplot(dataFrame_DR_70['X1'], shade=True, color="orange")
    ax = sns.kdeplot(dataFrame_DR_90['X1'], shade=True, color="green")
    plt.legend(['$80^\circ$', '$70^\circ$', '$90^\circ$'])
    plt.xlim([-2, 3])
    plt.xlabel("$X_1$")
    plt.show()
    #plt.savefig(method + '_features_pdf_' + str(angle) + 'deg.png')
    
    """
    p_00deg = plotDensity(dataFrame_DR_00['X1'][::1], dataFrame_DR_00['X2'][::1], 0)
    p_10deg = plotDensity(dataFrame_DR_10['X1'][::1], dataFrame_DR_10['X2'][::1], 10)
    p_20deg = plotDensity(dataFrame_DR_20['X1'][::1], dataFrame_DR_20['X2'][::1], 20)
    p_30deg = plotDensity(dataFrame_DR_30['X1'][::1], dataFrame_DR_30['X2'][::1], 30)
    p_40deg = plotDensity(dataFrame_DR_40['X1'][::1], dataFrame_DR_40['X2'][::1], 40)
    p_50deg = plotDensity(dataFrame_DR_50['X1'][::1], dataFrame_DR_50['X2'][::1], 50)
    p_60deg = plotDensity(dataFrame_DR_60['X1'][::1], dataFrame_DR_60['X2'][::1], 60)
    p_70deg = plotDensity(dataFrame_DR_70['X1'][::1], dataFrame_DR_70['X2'][::1], 70)
    p_80deg = plotDensity(dataFrame_DR_80['X1'][::1], dataFrame_DR_80['X2'][::1], 80)
    p_90deg = plotDensity(dataFrame_DR_90['X1'][::1], dataFrame_DR_90['X2'][::1], 90)
    """
    
    