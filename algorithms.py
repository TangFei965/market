import numpy as np
import scipy
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from fitsne import FItSNE

import pandas as pd

import scanpy as sc
import anndata


def normalize_to_TPM(array):
    
    """
    Parameters:
    -------------

    array: np.array()
        raw counts array

    Returns:
    -------------
    array: np.array()
        TPM array
    """

    #for row in range(0, array.shape[0]):
    #    row_sum = np.sum(array[row])
    #    array[row] = array[row] / row_sum * 1000000
    
    TPM = array / np.sum(array, axis=1).reshape(-1, 1) * 1e6

    return TPM

def select_informative_genes(TPM, genes, number_of_genes, labels = []):
    
    """
    Parameters:
    ---------------
    
    TPM: np.array()
        TPM of genes. cell_numbers by gene_numbers.
        
    genes: np.array()
        all genes to select from
        
    number_of_genes: int()
        number of genes to select
        
    Returns:
    ---------------
    
    informative_genes: np.array()
    informative_t_scores: np.array()
    informative_TPM: np.array()
    p_values: np.array()
    """
    
    # unsupervised gene selection
    if len(labels) == 0:
        labeled_TPM = TPM
    
    # supervised gene selection
    elif len(labels) != 0:
            
        groups = set(labels)
        group_TPMs = []

        for group in groups:
            group_indexes = [index for index, x in enumerate(labels) if x == group]
            group_TPM = TPM[group_indexes, :]
            group_TPM = np.mean(group_TPM, axis=0)
            group_TPMs.append(group_TPM)
            
        labeled_TPM = np.vstack(group_TPMs)
        # print('\n\n\n\n\n\n', labeled_TPM.shape, '\n\n\n\n')

    # calculate t_scores    
    TPM_1 = labeled_TPM + 1    
    log_E = np.log2(np.mean(TPM_1, axis=0))
    E_log = np.mean(np.log2(TPM_1), axis=0)
    t_scores = log_E - E_log
    
    # return values
    descending_indexes = t_scores.argsort()[::-1][:number_of_genes]
    informative_genes = genes[descending_indexes]
    informative_t_scores = t_scores[descending_indexes]
    informative_t_scores = np.around(informative_t_scores, decimals=1)

    informative_genes_indexes = []
    for i in informative_genes.tolist():
        informative_genes_indexes.append(genes.tolist().index(i))
    informative_TPM = TPM[:, informative_genes_indexes]

    p_vals = []
    if len(labels) != 0:
        label_sets = set(labels)
        for i in range(number_of_genes):
            expression_vector = informative_TPM[:,i]
            grouped_expressions = []
            for label in label_sets:
                grouped_expression = [expression_vector[i] for i in range(len(expression_vector)) if labels[i]==label]
                grouped_expressions.append(grouped_expression)
            p_vals.append(scipy.stats.f_oneway(*grouped_expressions)[1])
    p_vals = np.array(p_vals)

    return informative_genes, informative_t_scores, informative_TPM, p_vals

def tsne(X, method='umap'):
    
    """
    Parameters:
    ---------------

    X: numpy.array()
        log2(TPM + 1)
    method: str()
        'sklearn', 'FIt-SNE'

    Returns:
    ---------------
    
    tsne_1: list()
    tsne_2: list()
    """

    if method == 'sklearn':
        X_embedded = TSNE(n_components=2, metric='cosine', init='pca').fit_transform(X)

    # fit-sne seems to have the problem 'np.array() is not C - contiguous'
    elif method == 'FIt-SNE':
        # init
        pca = PCA(n_components=50)
        X_pca = pca.fit_transform(X)
        init = X_pca[:, :2]
        # make array C - contiguous
        X_pca = np.ascontiguousarray(X_pca)
        init = np.ascontiguousarray(init)
        # embedding
        X_embedded = FItSNE(X_pca, initialization=init)
    
    elif method == 'TSNEApprox': # from scanorama, the result seems to be very weird, do not suggest using.
        from scanorama.t_sne_approx import TSNEApprox
        # init
        pca = PCA(n_components=50)
        X_pca = pca.fit_transform(X)
        init = X_pca[:, 2]
        # embedding
        X_embedded = TSNEApprox(init='pca', metric='cosine').fit_transform(X_pca)
    
    elif method == 'umap':
        import umap
        reducer = umap.UMAP()
        X_embedded = reducer.fit_transform(X)

    tsne_1 = X_embedded[:,0].tolist()
    tsne_2 = X_embedded[:,1].tolist()

    return tsne_1, tsne_2


def calculate_unsupervised_clustering(TPM, genes, number_of_clusters, genes_for_clustering, method = 'louvain'):

    """
    Parameters:
    ---------------
    
    TPM: np.array()
        TPM of genes. cell_numbers by gene_numbers.
        
    genes: np.array()
        all genes to select from
        
    number_of_clusters: int()
        number of clusters in K-means

    genes_for_clustering: list()
        selected feature for clustring
        
    Returns:
    ---------------

    cell_types: list()
    """
    
    informative_genes_indexes = []
    for i in genes_for_clustering:
        informative_genes_indexes.append(genes.tolist().index(i))
    informative_TPM = TPM[:, informative_genes_indexes]

    if method == 'Kmeans':
        kmeans = KMeans(n_clusters=number_of_clusters).fit(informative_TPM)
        # kmeans = KMeans(n_clusters=number_of_clusters).fit(np.log2(informative_TPM))
        cell_types = kmeans.labels_.tolist()

    elif method == 'louvain':
        my_anndata = anndata.AnnData(informative_TPM)
        sc.tl.pca(my_anndata)
        sc.pp.neighbors(my_anndata)
        sc.tl.louvain(my_anndata)
        cell_types = my_anndata.obs['louvain'].tolist()

    elif method == '':
        pass

    def _descritize_cell_types(cell_types):
        for i in range(0, len(cell_types)):
            cell_types[i] = str(cell_types[i])
        return cell_types

    cell_types = _descritize_cell_types(cell_types)

    return cell_types
    
def convert_cross_species_genes(genes, conversion_reference):
    """
    Parameters:
    -----------
        genes: list()
        
        conversion_reference: dict()
    
    Returns:
    -----------
        converted_genes: list()
    """

    converted_genes = []
    for gene in genes:
        try:
            gene = conversion_reference[gene]
        except KeyError:
            gene = gene
        converted_genes.append(gene)
    
    return converted_genes

def integrate_datasets(datasets, genes, method='scanorama', merge=True):
    """
    Parameters:
    -----------
        datasets: list of scipy.sparse.csr_mstrix or np.array()

        genes: list of lists
    Returns:
    -----------
        integrated_dim_red: np.array() cell by dimensions
        corrected_expression_matrix: np.array() if merge, or lists of np.array() or scipy.sparse.csr_matrix
        corrected_genes: np.array()
    """
    if method == 'scanorama':
        import scanorama
        integrated_dim_red, corrected_expression_matrix, corrected_genes = \
            scanorama.correct(datasets, genes, dimred=100, return_dimred=True)
    
    if merge:
        corrected_expression_matrix = scipy.sparse.vstack(corrected_expression_matrix).toarray()
        integrated_dim_red = np.vstack(integrated_dim_red)
    
    return integrated_dim_red, corrected_expression_matrix, corrected_genes


# this should be put in db_query? 
def query_gene_expression_value(TPM, genes, user_input_genes, use_sum):
    """
    Parameters:
    ---------------
    
    TPM: np.array()
        TPM of genes. cell_numbers by gene_numbers.
        
    genes: np.array()
        all genes to select from
        
    user_input_genes: list()
        genes user input

    use_sum: Boolean
        whether or not use sum for gene expression value
        
    Returns:
    ---------------

    expression_vectors: nested list()
        in log2(TPM+1)
    """

    informative_genes_indexes = []
    for i in user_input_genes:
        informative_genes_indexes.append(genes.tolist().index(i))
    informative_TPM = TPM[:, informative_genes_indexes]

    expression_matrix = np.log2(informative_TPM + 1)

    if use_sum:
        expression_matrix = np.sum(expression_matrix, axis=1)
    
    expression_vectors = expression_matrix.T[0,:].tolist()

    return expression_vectors

#===================================================================

class ScibetClassifier:
    
    def __init__(self):
        self.reference_core = None
        self.reference_genes = None
        self.reference_cell_types = None
        
    def calculate_core(self, expr_matrix=None, genes=None, cell_types=None, 
                       select_genes=False, number_of_genes=1000, log_transformed=False):
        """
        Parameters:
        ----------
            expr_matrix: np.array() 
                    library size normalized
                    
            genes: np.array()
            
            cell_types: np.array()
            
            select_genes: bool()
                    if True, select informative genes before calculate scibet core,
            
            number_of_genes: int()
                    number of informative genes used for scibet core, if select_genes == True
                    
            log_transformed: bool()
                    final input for building scibet core should be log-transformed expression value.
                    if not log_transformed, perform log2(X+1) transformation
            
        Returns:
        ----------
            reference_core: np.array()
                log2(X) transformed probability
                
            reference_genes: np.array()
                
            reference_cell_types: np.array()
                labels
        """
        
        # select informative genes if select_genes == True:
        if select_genes:
            genes, _, expr_matrix, _, = \
                select_informative_genes(expr_matrix, genes, number_of_genes, labels=cell_types)
        
        # use log_transformed expression value
        if not log_transformed:
            expr_matrix = np.log2(expr_matrix + 1)
        
        # use pd.Dataframe for easier calculations
        df_expr_matrix = pd.DataFrame(expr_matrix)
        df_expr_matrix.index = cell_types
        df_expr_matrix.columns = genes
        
        # building df_prob
        df_prob = pd.DataFrame()
        cell_types_set = set(cell_types)
        for cell_type in cell_types_set:
            df_cell_type = df_expr_matrix[df_expr_matrix.index==cell_type]
            num_of_genes = df_cell_type.shape[1]
            expression_sum = np.sum(df_cell_type.to_numpy()) + num_of_genes
            
            cell_type_vector = []
            for gene in df_cell_type.columns:
                prob = (np.sum(df_cell_type[gene]) + 1) / expression_sum
                cell_type_vector.append(prob)    
            df_prob[cell_type] = cell_type_vector  
        df_prob.index = genes
        df_prob = df_prob.T
        
        df_prob_log = np.log2(df_prob)
        
        self.reference_core = df_prob_log.to_numpy()
        self.reference_genes = df_prob_log.columns.to_numpy()
        self.reference_cell_types = df_prob_log.index.to_numpy()

        return self.reference_core, self.reference_genes, self.reference_cell_types
    
    def predict_cell_type(self, expr_matrix=None, genes=None, log_transformed=False):
        """
        Parameters:
        -----------
            expr_matrix: np.array()
                     library size normalized 
        
            genes: np.array()
            
            log_transformed: bool()
                    final input for building scibet core should be log-transformed expression value.
                    if not log_transformed, perform log2(X+1) transformation
        
        Returns:
        -----------
            predicted_cell_types: np.array()

            df_prob: pd.DataFrame()
                    query cells by reference cells probability dataframe
            
        """
        
        # use log_transformed expression value
        if not log_transformed:
            expr_matrix = np.log2(expr_matrix + 1)
        
        df_expr_matrix = pd.DataFrame(expr_matrix)
        df_expr_matrix.columns = genes
        df_expr_matrix_informative_genes = df_expr_matrix.reindex(columns=self.reference_genes, fill_value=0)
        
        # no need to transform expression level of test matrix to prob, as the final result of argmax(y) should be the same
        df_pred = pd.DataFrame(np.dot(df_expr_matrix_informative_genes, self.reference_core.T), columns=self.reference_cell_types)
        
        predicted_cell_types = df_pred.T.idxmax(0).to_numpy()
        self.predicted_cell_types = predicted_cell_types
        
        prob_array = df_pred.to_numpy()
        prob_array = 2 ** (prob_array - np.max(prob_array, axis=1).reshape(-1,1))
        
        df_prob = pd.DataFrame(prob_array / np.sum(prob_array, axis=1).reshape(-1, 1), index=df_pred.index, columns=df_pred.columns)
        
        return self.predicted_cell_types, df_prob 
    
    def generate_projection(self, reference_cor_1, reference_cor_2, reference_cell_types):
        """
        Parameters:
        -----------
            reference_cor_1: np.array()
            
            reference_cor_2: np.array()
            
            reference_cell_types: np.array()
        
        Returns:
        ----------
            projection_cor_1: list()
            
            projection_cor_2: list()
        """
        
        projection_cor_1 = []
        projection_cor_2 = []
        
        import random
        for cell_type in self.predicted_cell_types:
            indexes = [i for i, x in enumerate(reference_cell_types) if reference_cell_types[i] == cell_type]
            random_index = random.choice(indexes)
            projection_cor_1.append(reference_cor_1[random_index]) #+ float(10*np.random.rand(1) * 0.2 - 0.1))
            projection_cor_2.append(reference_cor_2[random_index]) # + float(10*np.random.rand(1) * 0.2 - 0.1))
        
        return projection_cor_1, projection_cor_2
    
    def generate_heatmap(self):
        """
        draw a clustered heatmap to show marker genes of reference cell types

        Parameters:
        -----------

        Returns:
        -----------
            heatmap_cell_types: np.array()

            heatmap_genes: np.array()

            heatmap_matrix: np.array()
        """
        print('before seaborn')
        import seaborn as sns
        print('seaborn')
        print(self.df_TPM_averaged.shape)

        df_heatmap = np.log2(self.df_TPM_averaged + 1)
        print(df_heatmap.shape)
        my_cluster_map = sns.clustermap(df_heatmap)
        print('finished ploty')

        row_index = my_cluster_map.dendrogram_row.reordered_ind
        column_index = my_cluster_map.dendrogram_col.reordered_ind
        print(row_index)
        print(column_index)

        df_cluster_map = df_heatmap.iloc[row_index, column_index]
        print(df_cluster_map.shape)

        heatmap_cell_types = df_cluster_map.index.to_numpy()
        heatmap_genes = df_cluster_map.columns.to_numpy()
        heatmap_matrix = df_cluster_map.to_numpy()

        return heatmap_cell_types, heatmap_genes, heatmap_matrix

