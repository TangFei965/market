import numpy as np
import pandas as pd


class toTPM():

    def __init__(self, data: np.array=None, genes: list=None, species=None):
        """
        Parameters:
        -----------
        genes: np.array()
            genes names
            n genes
    
        """
        self.data = data
        self.genes = genes
        self.species = species

        

    def read_length_from_db(self):
        """
        load the basic length from  the annotations
        """

        if self.species=='human':
            gene_length_list = pd.read_csv('/home/biodb/data/dataset_collection/resources/gene_references/gene_length/human_unique.tsv', sep ='\t', index_col=1)

        if self.species=='mouse':
            gene_length_list = pd.read_csv('/home/biodb/data/dataset_collection/resources/gene_references/gene_length/mouse_unique.tsv', sep ='\t', index_col=1)

        gene_length_dic = {}
        for i in gene_length_list.index:
            gene_length_dic[i]=gene_length_list.loc[i,'Gene_length']

        return gene_length_dic

    
    def get_data_gene_length(self):
        """
        get the target length for aimed genes
        """

        #count_list=[]
        #for i in genes:
        #    count_list.append(gene_length_dic[i])
        length_list = [gene_length_dic[i] for i in self.genes]

        return count_list


    
    def count_TPM(self):
        """
        count TPM for data
        """

        new_data = np.zeros(self.data.shape)
        for j in range(self.data.shape[1]): 
            sequence_depth = self.data.sum(axis=0)[j]
            for i in range(len(self.genes)):
                new_data[i][j] = ( i * pow(10,6) ) / ( length_list[i] * sequence_depth )

        return new_data







