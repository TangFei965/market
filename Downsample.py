import numpy as np
import pandas as pd

import seaborn as sea
import matplotlib.pyplot as plt
from collections import Counter

import scanpy as sc
import anndata

class Downsample():
    
    def __init__(self, adata=None):
        """
        Parameters:
        -----------
        adata: two-dimensional np.array()
            t-sne cor

        -------------
        array: np.array()
        
        """
        self.adata = adata

    def downsample(self):

        df_tsne_cor=pd.DataFrame(adata, columns=['TSNE-1','TSNE-2'])
        raw_cor = df_tsne_cor.to_numpy()
        floor_cor = np.floor(raw_cor)

        grids = np.array([str(floor_cor[i][0]) + ' * ' + str(floor_cor[i][1]) for i in range(0, floor_cor.shape[0])])
        df_tsne_cor['grid'] = grids
        counts = Counter(grids)

        sampled_tsne_1 = np.array([]).reshape(-1, 1)
        number_in_grid = []
        for grid in set(grids):
            # print('%sth itetration' %i)
            indexes = np.where(grids == grid)
            grid_cor = raw_cor[indexes]

            # sampling number please see here.
            sampling_number = np.floor(np.log(counts[grid] + 1)).astype(int)
            _tsne_1 = np.random.choice(grid_cor[:, 0], sampling_number, replace=False)
            number_in_grid.append(_tsne_1.shape[0])
            sampled_tsne_1 = np.append(sampled_tsne_1, _tsne_1.reshape(-1,1), axis=0)
        
        sampled_tsne_1=sampled_tsne_1.reshape(len(sampled_tsne_1))
        df_tsne_cor_sampled = df_tsne_cor[df_tsne_cor['TSNE-1'].isin(sampled_tsne_1)]
        sampled_corrs = df_tsne_cor_sampled.to_numpy()

        return sampled_corrs

    def downsampleGraph(self, sampled_corrs):
        plt.figure(figsize=(15,15))
        plt.scatter(sampled_corrs[:,0], sampled_corrs[:,1], s=5)
        
        return sampled_corrs

    def filtercell(adata):
        # Purify cell population
        filter_data=adata[df_tsne_cor_sampled.index]

        sc.pp.neighbors(filter_data, n_neighbors=20, n_pcs=40)
        sc.tl.umap(filter_data)
        sc.tl.louvain(filter_data)

        raw_adata=adata
        filter_data.obs['raw_louvain']=raw_adata[filter_data.obs.index].obs['louvain']

        ad=filter_data.obs

        a={}
        for i in range(23):
            a[i]=Counter(ad[ad['louvain']=='%s' %i].loc[:,'raw_louvain'])

        b={}
        for i in range(32):
            b[i]=Counter(ad[ad['raw_louvain']=='%s' %i].loc[:,'louvain'])
        
        cell_index=[]
        for i in range(len(a)):
            for j in a[i]:
                p1 = a[i][j] / sum(a[i].values())
                if p1 > 0.8:
                    cell_index.append(ad[(ad['raw_louvain']==j)&(ad['louvain']==i)].index)
                else:
                    p2 = b[j][i] / sum(b[j].values())
                    if p2 > 0.9:
                        cell_index.append(ad[(ad['raw_louvain']==j)&(ad['louvain']==i)].index)
        
        h=[]
        for i in range(32):
            h += cell_index[i].values.tolist()

        adata_filter=adata[h]

        return adata_filter

    def clusterCorrelation(adata_filter):
        # calculate the Correlation between cell clusters
        ds_TM_filter=[]
        test={}
        for i in range(23):
            ds_TM_filter += [adata_filter[ad_filter['louvain']==i].X.todense().mean(axis=0)]
            test[i]=set(adata_filter[adata_filter.obs['louvain']==i].loc[:,'raw_louvain'])

        raw_ds_TM=[]
        for i in range(len(ds_TM_filter)):
            for j in test[i]:
                raw_ds_TM += [raw_adata[raw_adata.obs['louvain']=='%s' %j].X.todense().mean(axis=0)]

        ds_TM_filter=pd.DataFrame(np.array(ds_TM_filter).reshape(23,1000).T)
        raw_ds_TM_filter=pd.DataFrame(np.array(raw_ds_TM).reshape(32,1000).T)

        coef, p = spearmanr(ds_TM_filter, raw_ds_TM_filter)

        sea.boxplot(coef)

        return coef, p

                

        















        







    

    
   




