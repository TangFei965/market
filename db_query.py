from pymongo import MongoClient
import numpy as np

import re
##################################

# this should be further uncoupled into building block functions

def get_all_public_datasets():
    """
    Parameters:
    -----------
    
    Returns:
        
        public_datasets: list()
    """
    client = MongoClient()
    dbs = client.list_database_names()
    public_datasets = []
    for db_name in dbs:
        if db_name.startswith('GSE'):
            public_datasets.append(db_name)
    
    return public_datasets

def get_tsne(GSE_number):
    """
    Parameters:
    -----------
        GSE_number: str()
    
    Returns:
    -----------
        tsne_1: list()
        tsne_2: list()
        cell_types: list()
    """
    client = MongoClient()
    db = client[GSE_number]

    tsne_1 = db['pre_calculated_tsne'].find({"field": "tsne_1"}).next()['content']
    tsne_2 = db['pre_calculated_tsne'].find({"field": "tsne_2"}).next()['content']
    cell_types = db['pre_calculated_tsne'].find({"field": "cell_type"}).next()['content']

    return tsne_1, tsne_2, cell_types

def get_tsne_center(GSE_number):
    """
    Parameters:
    -----------
        GSE_number: str()
    
    Returns:
    -----------
        tsne_1_center: list()
        tsne_2_center: list()
        cell_types_center: list()
    """

    from statistics import mean, median

    tsne_1, tsne_2, cell_types = get_tsne(GSE_number)
    tsne_1_center = []
    tsne_2_center = []
    cell_types_center = []
    for cell_type in set(cell_types):
        cell_type_indexes = [index for index, x in enumerate(cell_types) if x==cell_type]
        cell_type_tsne_1 = [tsne_1[i] for i in cell_type_indexes]
        cell_type_tsne_2 = [tsne_2[i] for i in cell_type_indexes]
        tsne_1_center.append(median(cell_type_tsne_1))
        tsne_2_center.append(median(cell_type_tsne_2))
        cell_types_center.append(cell_type)
    
    return tsne_1, tsne_2, cell_types, tsne_1_center, tsne_2_center, cell_types_center

def get_genes(GSE_number, type='symbol'):
    """
    Parameters:
        GSE_number: str()
        
        type: str()
            gene type, either 'symbol', 'ensemble', 'GeneID'
    Returns:
        genes: list()
    """
    
    client = MongoClient()
    db = client[GSE_number]
    collection = db['gene']
    genes = collection.find_one({"type": type})['content']

    return genes

##################################
##################################


"""
use the command line command to import json

mongoimport -d 00000 --collection matrix --drop --file 2000_cells.json 

"""

"""

design schema 

db designs

db: doi:00000

    collections:

        info:   this is a collection(table)
            for paper info

        matrix: this is another collection(table)
            for the matrix

            {"genes": ['a', 'b', 'c']},

            x n elements      
                { 
                    "_id" : ObjectId("542c2b97bac0595474108b48"),  system generated
                    
                    "index": 0,

                    "matrix": {

                        "labels":{
                            "cell_type": ...
                        }

                        "data":[] list. need to be changed into sparse matrix, but csr is not json serializable...
                    }
                
                }


"""



def query_index_of_all_datasets():
    """
    response = {
        'index_of_all_datasets': [
            {
                'title': 'RNA Sequencing of Single Human Islet Cells Reveals Type 2 Diabetes Genes',
                'GSE_number': 'GSE81608'
            },
            {
            },
        ]
    }
    """

    client = MongoClient()
    dbs = client.list_database_names()
    public_datasets = []
    for db_name in dbs:
        if db_name.startswith('GSE'):
            public_datasets.append(db_name)

    index_of_all_datasets = []
    errored_datasets = []
    for public_dataset in public_datasets:
        try:
            db = client[public_dataset]
            title = db['info'].find({"field": "title"}).next()['content']
            abstract = db['info'].find({"field": 'abstract'}).next()['content']

            try: 
                figure_url = db['info'].find({"field": "figure_url"}).next()['content']
            except:
                try:
                    figure_url = db['info'].find({"field": "_figure_url"}).next()['content']
                except:
                    figure_url = ''
            
            index_of_all_datasets.append({
                'title': title,
                'abstract': abstract,
                'GSE_number': public_dataset,
                'figure_url': figure_url,
            })
        except Exception as e:
            errored_datasets.append({"dataset_name": public_dataset, "error_message": e})
            
    print('errored_datasets', errored_datasets)
            
    number_of_datasets = len(index_of_all_datasets)

    response = {
        "number_of_datasets": number_of_datasets,
        'index_of_all_datasets': index_of_all_datasets,
    }   

    return response

def get_taxonomy_id(GSE_number):
    """
    Parameters:
    -----------
        GSE_number: str()
        
    Returns:
    -----------
        taxonomy_id: str()
    """
    
    client = MongoClient()
    db = client[GSE_number]
    collection = db['info']
    
    taxonomy_id = collection.find_one({"field": "taxonomy_id"})['content']
    
    return taxonomy_id

def get_gene_conversion_reference(taxonomy_id = '10090', conversion_type = 'ensemble_to_symbol'):
    """
    Parameters:
    -----------
        taxonomy_id: str()
            species id
        conversion_type: str()
            gene type conversion type
    Returns:
    -----------
        conversion_reference: dict()
    """
    
    client = MongoClient()
    db = client['GENE']
    collection = db['tax_id_' + taxonomy_id]
    conversion_reference = collection.find({"field": conversion_type}).next()['content']

    return conversion_reference

def get_species_gene_conversion_reference(input_species = '10090', output_species = '9606'):
    """
    Parameters:
    -----------
        input_species: str()
            default = mouse
        output_species: str()
            default = human
    Returns:
    ------------
        conversion_reference: dict()
    """

    client = MongoClient()
    db = client['GENE']
    collection = db['conversion']
    conversion_reference = collection.find_one({"type": input_species + '_to_' + output_species})['content']

    return conversion_reference

def get_all_genes_of_species(taxonomy_id='10090', type='symbol'):
    """
    Parameters:
    -----------
        taxonomy_id: str()
            species id
        type: str()
            either 'symbol' or 'ensemble'
    Returns:
    -----------
        species_genes: list()
    """

    client = MongoClient()
    db = client['GENE']
    collection = db['tax_id_' + taxonomy_id]
    species_genes = collection.find_one({"field": "all_" + type})['content']

    return species_genes

def query_individual_dataset(doi):
    """
    returns the basic info of an individual dataset 
    """

    client = MongoClient()
    db = client[doi]
    print(doi)
    
    # info collection
    title = db['info'].find({"field": "title"}).next()['content']
    doi = db['info'].find({"field": "doi"}).next()['content']
    abstract = db['info'].find({"field": 'abstract'}).next()['content']
    
    try:
        doi_link = db['info'].find({"field": "doi_link"}).next()['content']
    except:
        doi_link = ''

    try: 
        figure_url = db['info'].find({"field": "figure_url"}).next()['content']
    except:
        try:
            figure_url = db['info'].find({"field": "_figure_url"}).next()['content']
        except:
            figure_url = ''

    # pre_calculated_tsne_collection
    tsne_1 = db['pre_calculated_tsne'].find({"field": "tsne_1"}).next()['content']
    tsne_2 = db['pre_calculated_tsne'].find({"field": "tsne_2"}).next()['content']
    cell_type = db['pre_calculated_tsne'].find({"field": "cell_type"}).next()['content']

    # marker_genes
    collection = db['marker_genes']
    cursor = collection.find()
    markers = []
    while True:
        try:
            document = cursor.next()
            marker = {
                "cell_type": document['cell_type'], 
                "marker_genes": document["marker_genes"], 
                "marker_t_scores": document["marker_t_scores"]
            }
            markers.append(marker)
        except StopIteration:
            break 

    response = {

        "info":{
            "title": title,
            "abstract": abstract,
            "doi": doi,
            "doi_link": doi_link,
            "figure_url": figure_url,
        },

        "pre_calculated_tsne":{
            "tsne_1": tsne_1,
            "tsne_2": tsne_2,
            "cell_type":cell_type
        },

        "markers": markers

    }
    
    return response

def query_gene_expression_individual_dataset(GSE_number='', gene='', indexes=[]):
    """
    Paramerters:
    ------------
        GSE_number: str()
            dataset to query from
        gene: str()
            gene to query
        indexes: list()
            cells of interests, default to empty list, meaning all cells
    
    Returns:
    --------
        not_found: Boolean()
            True if not found,
        expression_vector: list()
            gene expression vector
    """

    client = MongoClient()
    db = client[GSE_number]
    collection = db['transposed_matrix']

    # original code without support of gene conversion
    #gene_list = db['matrix'].find_one()['genes'] # should be changed using .var collection in the future
    #if gene in gene_list:
    #    not_found = False
    #    print(collection.find({"gene": gene}), type(collection.find({"gene": gene})))
    #    expression_vector = collection.find({"gene": gene}).next()['data']
    #    expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
    #else:
    #    not_found = True
    #    expression_vector = []

    # new code to support gene conversion
    # this should be further revised after the gene format is finally dertermined.
    gene_name_collection = db['gene']
    not_found = False
    taxonomy_id = get_taxonomy_id(GSE_number)
    
    if gene.startswith('ENS'):
        print(gene)
        try:
            expression_vector = collection.find({"gene": gene}).next()['data']
            expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
        except:
            try:
                conversion_reference = get_gene_conversion_reference(taxonomy_id=taxonomy_id, \
                    conversion_type='ensemble_to_symbol')
                gene = conversion_reference[gene]
                print('converted_gene', gene)
                expression_vector = collection.find({"gene": gene}).next()['data']
                expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
            except:
                print('reached here')
                not_found = True
                expression_vector = []                          
    elif gene[0].isdigit():
        try:
            conversion_reference = get_gene_conversion_reference(taxonomy_id=taxonomy_id, \
                conversion_type='GeneID_to_symbol')
            converted_gene = conversion_reference[gene]
            print('converted_gene', converted_gene)
            expression_vector = collection.find({"gene": converted_gene}).next()['data']
            expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
        except:
            try:
                conversion_reference = get_gene_conversion_reference(taxonomy_id=taxonomy_id, \
                    conversion_type='GeneID_to_ensemble')
                converted_gene = conversion_reference[gene]
                print('converted_gene', converted_gene)
                expression_vector = collection.find({"gene": converted_gene}).next()['data']
                expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
            except:
                print('reached here')
                not_found = True
                expression_vector = []
    else: # input is gene symbol
        print(gene)
        try:
            expression_vector = collection.find({"gene": gene}).next()['data']
            expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
        except:
            try: # case insensitive gene symbol search
                print('case insensitive gene symbol search worked')
                genes_in_gene_symbol = get_genes(GSE_number)
                candidate_genes = [i for i in genes_in_gene_symbol if re.search(gene, i, re.IGNORECASE)]
                print(candidate_genes)
                case_converted_gene = candidate_genes[0]
                print(case_converted_gene)
                expression_vector = collection.find({"gene": case_converted_gene}).next()['data']
                expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
            except:
                try:
                    conversion_reference = get_gene_conversion_reference(taxonomy_id=taxonomy_id, \
                        conversion_type='symbol_to_ensemble')
                    gene = conversion_reference[gene]
                    expression_vector = collection.find({"gene": gene}).next()['data']
                    expression_vector = np.log2(np.array(expression_vector) + 1).tolist()
                except:
                    print('reached here')
                    not_found = True
                    expression_vector = []       
    return not_found, expression_vector

def query_gene_expression_across_all_datasets(gene_of_interest, page):
    """
    Parameters:
        gene_of_interest: str()
    --------------

    Returns:
        response = {

            'number_of_datasets': '',
            
            'datasets':[
                    {
                    'title': '',
                    'GSE': '',
                    'expression_vector': '',
                    'cell_type': [],
                }            
            ],

            'end_of_datasets': '',  # if end of datasets, 'end_of_datasets': 'end_of_datasets'

        }
    -------------
    """

    client = MongoClient()

    # first show all public datasets
    public_datasets = get_all_public_datasets()

    # then get all datasets which has the gene of interest
    datasets_with_gene_of_interest = []
    for dataset in public_datasets:
        db = client[dataset]
        collection = db['matrix']
        all_genes = collection.find_one()['genes']
        if gene_of_interest in all_genes:
            datasets_with_gene_of_interest.append(dataset)
    
    number_of_datasets = len(datasets_with_gene_of_interest)
    
    # the front end will request 'a page of dataset, currenlt 5 dataset each page' each time
    # now we determine which page we should be on per request
    current_page_datasets = datasets_with_gene_of_interest[0+5*(page-1): 4+5*(page-1)]
    end_of_all_datasets = ''
    if len(current_page_datasets) < 5:
        end_of_all_datasets = 'end_of_all_datasets'
    
    # actually traversing the dataset and querying expression value in log2(TPM+1)
    _datasets = []
    for dataset in current_page_datasets:
        
        db = client[dataset]
        
        title = db['info'].find({"field": "title"}).next()['content']
        GSE = dataset
        
        ########========================== original code, deprecated  ==========================########
        #collection = db['matrix']
        #cursor = collection.find()
        #all_genes = cursor.next()['genes']
        #gene_index = all_genes.index(gene_of_interest)
        
        #expression_vector = []
        #cell_types = []
        #while True:
        #    try:
        #        cell = cursor.next()
        #        cell_type = cell['matrix']['labels']['cell_type']
        #        cell_types.append(cell_type)
        #        cell_expression_vector = np.array(cell['matrix']['data'])
        #        cell_expression_vector = cell_expression_vector / np.sum(cell_expression_vector) * 1000000
        #        expression_value = np.log2(cell_expression_vector.tolist()[gene_index] + 1)
        #        expression_vector.append(np.around(expression_value, 2))
        #        
        #    except StopIteration:
        #        break

        ##### new version code
        _, expression_vector = query_gene_expression_individual_dataset(GSE_number=GSE, gene=gene_of_interest, indexes=[])
        _, _, cell_types = get_tsne(GSE)

        _dataset = {
            'title': title,
            'GSE': GSE,
            'expression_vector': expression_vector,
            'cell_type': cell_types,
        }   
        
        _datasets.append(_dataset)

    response = {
        'number_of_datasets': number_of_datasets,
        'datasets': _datasets,
        'end_of_datasets': end_of_all_datasets,
        'page': page
    }

    return response

def query_cells_of_given_indexes(doi, indexes):
    """
    receives a list of cell indexes (starting from 0), and get the info of these cells
    
    Parameter:
    ----------
    doi: 
    
    indexes: list()
    
    Returns:
    currently a np.array(), in the future rewrite with decorators
    ----------
    
    """
    
    client = MongoClient()
    db = client[doi]
    collection = db['matrix']

    rows = []

    for index in indexes:
        cursor = collection.find({"index": index})
        
        row = cursor.next()['matrix']['data']
        rows.append(row)
    
    genes = np.array(get_genes(doi))
    
    TPM = np.vstack(rows)
    
    return TPM, genes   

def query_scibet_core(dataset_id):
    """
    Parameters:
    -----------
        dataset_id: str()
            the dataset id

    Returns:
    -----------
        reference_core: np.array()
            log2(X + 1), X being normalized probability

        reference_genes: np.array()

        reference_cell_types: np.array()
        
    """

    client = MongoClient()
    db = client[dataset_id]
    collection = db['scibet_core']
    
    documents = collection.find({})
    reference_cell_types = []
    probability_vectors = []
    reference_genes = []
    for document in documents:
        try:
            reference_cell_types.append(document['cell_type'])
            probability_vectors.append(document['prob_vector'])
        except:
            reference_genes = document['reference_genes']
    
    reference_core = np.array(probability_vectors)
    reference_cell_types = np.array(reference_cell_types)
    reference_genes = np.array(reference_genes)
    
    return reference_core, reference_cell_types, reference_genes
            
    