

def downsample(adata):

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

        return df_tsne_cor_sampled

    def filtercell(filter_data, adata):
        # Purify cell population
      
        sc.pp.neighbors(filter_data, n_neighbors=20, n_pcs=40)
        sc.tl.umap(filter_data)
        sc.tl.louvain(filter_data)

        raw_adata=adata
        filter_data.obs['raw_louvain']=raw_adata[filter_data.obs.index].obs['louvain']

        ad=filter_data.obs

        a={}
        for i in range(len(set(filter_data.obs['louvain']))):
            a[i]=Counter(ad[ad['louvain']=='%s' %i].loc[:,'raw_louvain'])

        b={}
        for i in range(len(set(raw_adata.obs['louvain']))):
            b[i]=Counter(ad[ad['raw_louvain']=='%s' %i].loc[:,'louvain'])
        
        cell_index=[]
        for i in range(len(a)):
            for j in a[i]:
                p1 = a[i][j] / sum(a[i].values())
                if p1 > 0.7:
                    cell_index.append(ad[(ad['raw_louvain']==j)&(ad['louvain']=='%s' %i)].index)
                else:
                    p2 = b[int(j)][i] / sum(b[int(j)].values())
                    if p2 > 0.7:
                        cell_index.append(ad[(ad['raw_louvain']==j)&(ad['louvain']=='%s' %i)].index)
        
        h=[]
        for i in range(len(cell_index)):
            h += cell_index[i].values.tolist()

        adata_filter=filter_data[h]
        return adata_filter

    def clusterCorrelation2(adata_filter,raw_adata):
        ds_corr=[]
        for i in set(adata_filter.obs['louvain']):
            ds_filter=[]
            test={}
            raw_ds=[]
            ds_filter += [adata_filter[adata_filter.obs['louvain']=='%s' %i].X.todense().mean(axis=0)]
            test[i]=set(adata_filter[adata_filter.obs['louvain']=='%s' %i].obs.loc[:,'raw_louvain'])
            for j in test[i]:
                raw_ds += [raw_adata[raw_adata.obs['louvain']==j].X.todense().mean(axis=0)]
                ds_filter=pd.DataFrame({"ds_filter": np.array(ds_filter).reshape(np.array(ds_filter).shape[2]), 
                                        "raw_ds": np.array(raw_ds).reshape(np.array(raw_ds).shape[2])})
                #calculate the corelation
                ds_corr.append(np.triu(ds_filter.corr(method="spearman"))[0][1])
        return ds_corr
    
