import scanpy as sc
import os
cwd = os.getcwd()
os.chdir('/home/biodb/data/Abio/Abio_test')
from market.load_db import db_loader
from market.algorithms import *
from market.db_query import *
os.chdir(cwd)
import pandas as pd
import numpy as np

#load data
GSE85241=pd.read_csv("/data/biodb/personal_codes/tf/GSE81076_D2_3_7_10_17.txt",sep='\t',index_col=0)
GSE81076=pd.read_csv("/data/biodb/personal_codes/tf/GSE85241_cellsystems_dataset_4donors_updated.csv",sep='\t', index_col=0)

#intergrate with scanorama
a = GSE85241.values.T
b = GSE81076.values.T
datasets=[a,b]
agene = GSE85241.index.tolist()
bgene = GSE81076.index.tolist()
genes = [agene,bgene]

pancreas_intergrate=integrate_datasets(datasets,genes)


NAMESPACE = 'pancreas'
data_names = [
    'GSE85241',
    'GSE81076'
]

if __name__ == '__main__':
    datasets, genes_list, n_cells = load_names(data_names)

label





