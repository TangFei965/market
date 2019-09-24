from django.http import HttpResponse
from django.http import JsonResponse

from django.views.decorators.csrf import csrf_exempt

import os
import json
import re
from pprint import pprint

from .db_query import *
from .algorithms import *
from .load_db import *

import pandas as pd
import numpy as np


#####################################################################################################################################


#  functions for specific requests


#####################################################################################################################################

def search_gene_across_all_datasets(order):
    """

     dataToServer= {

        "method": { 
            "type": "search_gene_across_all_datasets",
            "parameters": { 
                "gene_of_interest": 'erbb2', // str()
                "page": 1,  // int()
            }
        },
    
    },

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

    """
    gene_of_interest = order["method"]["parameters"]["gene_of_interest"]
    page = order["method"]["parameters"]["page"]
    response = query_gene_expression_across_all_datasets(gene_of_interest, page)

    return response

def merge_expression_matrix(order):
    """
    dataToServer= {

        // type of analysis
        "method": { 
            "type": "plot_tsne",
            "parameters": { "number_of_genes" : 500}      
        },
    
        "basket":{
            "basketID": [], use this to check if len() == 0
            "doi": ['00000', '00000', '00000', '00000'],
            "indexes": [0,1,2,3]  // indexes are intrisically self-contained in arrays.
        },
    
    },

    Returns:
    ---------
        dataset_composition: str()
            choose from list: ["single_dataset", "cross_datasets", "error"]

        error_message: str()
        
        candidate_TPM: np.array()

        candidate_genes: np.array()
    """
    
    """
    # test code
    basketID = []
    # basketID = ['a', 'a', 'a', 'a', 'a', 'a','a', 'a', 'a','a', 'a', 'a','a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b',  'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'a','a', 'a', 'a','a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b',  'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c']
    doi = '00000'
    indexes = list(range(100))
    number_of_genes = 200
    """
    
    doi = order['basket']['doi']
    indexes = order['basket']['indexes']

    datasets = list(set(doi))
    species = set([get_taxonomy_id(i) for i in datasets])
    print('\ndatasets:', datasets)
    print('\nspecies:', species)
    
    if len(datasets) == 1:
        dataset_composition = "single_dataset"
        error_message = ''
        doi = datasets[0]
        integrated_dim_red = None
        candidate_TPM, candidate_genes = query_cells_of_given_indexes(doi, indexes)
    
    elif len(datasets) > 1: #species == ['9606', '10090']:
        
        if species == set(['9606', '10090']) or species == set(['9606']) or species == set(['10090']):
            dataset_composition = "cross_datasets"
            error_message = ''

            # query each block and returns lists for scanorama integration
            TPM_lists = []
            genes_lists = []
            current_doi = doi[0]
            current_doi_indexes = []
            for i in range(len(doi)):
                _doi = doi[i]
                if _doi == current_doi and i != len(doi) - 1:
                    current_doi_indexes.append(indexes[i])
                elif _doi != current_doi or i == len(doi) - 1: # query each block and then append to the lists for scanorama
                    if i == len(doi) - 1:
                        current_doi_indexes.append(indexes[i])
                    current_TPM, current_genes, = query_cells_of_given_indexes(current_doi, current_doi_indexes)    
                    current_genes, _, current_TPM, _ = \
                        select_informative_genes(current_TPM, current_genes, 2000)
                    TPM_lists.append(current_TPM)
                    current_genes = current_genes.tolist()
                    if species == set(['9606', '10090']) and get_taxonomy_id(current_doi) == '10090':
                        conversion_reference = get_species_gene_conversion_reference() 
                        current_genes = convert_cross_species_genes(current_genes, conversion_reference)
                    genes_lists.append(current_genes)
                    
                    print('\t'*3, current_doi)
                    current_doi = _doi # start next block
                    current_doi_indexes = [indexes[i]]
                    
            # first try if scanorama works
            # get corrected expression
            integrated_dim_red, corrected_TPM, corrected_genes = integrate_datasets(TPM_lists, genes_lists)
            candidate_TPM, candidate_genes = corrected_TPM, corrected_genes
            
            # scanorama results are mostly between 0 and 1 and contains negative entries with very small absolute values \ 
            # we substract the minimal value (so that all counts are positive), normalize to TPM \
            # and treat this as a normal TPM
            
            np.seterr(all='raise')
            try:
                candidate_TPM = normalize_to_TPM(candidate_TPM - np.min(candidate_TPM, axis=1).reshape(-1, 1))
            except:
                print("scanorama fails, this is very likely that the datasets share no common cells")
                new_TPM_lists = []
                for i in range(len(TPM_lists)):
                    current_TPM = TPM_lists[i]
                    current_genes = genes_lists[i]
                    current_indexes = [i for i in range(len(current_genes)) if current_genes[i] in corrected_genes]
                    print(current_indexes)
                    current_TPM = current_TPM[:, current_indexes]
                    new_TPM_lists.append(current_TPM)
                candidate_TPM = normalize_to_TPM(np.vstack(new_TPM_lists))
    
        else:
            dataset_composition = 'error'
            error_message = 'Sorry, currenty we only support cross species analysis on mouse and human data.'
            integrated_dim_red, candidate_TPM, candidate_genes = None, None, None
            

    return dataset_composition, error_message, integrated_dim_red, candidate_TPM, candidate_genes

def plot_tsne(order, current_user):

    """
    dataToServer= {

        // type of analysis
        "method": { 
            "type": "plot_tsne",
            "parameters": { "number_of_genes" : 500}      
        },
    
        "basket":{
            "basketID": [], use this to check if len() == 0
            "doi": ['00000', '00000', '00000', '00000'],
            "indexes": [0,1,2,3]  // indexes are intrisically self-contained in arrays.
        },
    
    },
    """

    # real code
    basketID = order['basket']['basketID']
    number_of_genes = order['method']['parameters']['number_of_genes']
    #doi = set(order['basket']['doi']).pop()  # currently no cross datasets requests, use set() for simplicity

    #indexes = order['basket']['indexes']
    
    
    # candidate_TPM, candidate_genes = query_cells_of_given_indexes(doi, indexes)
    dataset_composition, error_message, integrated_dim_red, candidate_TPM, candidate_genes = \
         merge_expression_matrix(order)

    if dataset_composition == 'single_dataset':
        error = False
        informative_genes, informative_t_scores, informative_TPM, p_vals = \
                select_informative_genes(candidate_TPM, candidate_genes, number_of_genes, labels=basketID)
        # calculate tsne
        tsne_1, tsne_2 = tsne(np.log2(informative_TPM + 1))
    
    elif dataset_composition == 'cross_datasets':
        if candidate_genes.shape[0] < number_of_genes:
            number_of_genes = candidate_genes.shape[0]
        error = False
        informative_genes, informative_t_scores, informative_TPM, p_vals = \
                select_informative_genes(candidate_TPM, candidate_genes, number_of_genes, labels=basketID)
        tsne_1, tsne_2 = integrated_dim_red[:, 0].tolist(), integrated_dim_red[:, 1].tolist()

        # current_user.merged_matrix = candidate_TPM save for future use
    
    elif dataset_composition == 'error':
        error = True
        informative_genes, informative_t_scores, tsne_1, tsne_2 = np.array([]), np.array([]), [], []


    dataToClient = {
    
        "genes": [
                informative_genes.tolist(), # genes in ascending order
                informative_t_scores.tolist()    # t scores in ascending order
            ],     

        "matrix": {      
            "labels":{
                "tsne_1": tsne_1,
                "tsne_2": tsne_2, 
            }
        },

        "basketID": basketID,

        "message": {
            "error": error,
            "error_message": error_message,
        },

    }

    return dataToClient

def unsupervised_clustering(order, current_user):

    """
    dataToServer= {

        "method": { 
            "type": "unsupervised_clustering",
            "parameters": { 
                "number_of_clusters": 10,
                "genes_for_clustering": ['a', 'b', 'c', 'd']   // same genes used as in the TSNE graph
            }      
        },
    
        "basket":{
            "doi": ['00000', '00000', '00000', '00000'],
            "indexes": [0,1,2,3]
        },
    
    },
    """

    """
    # test code
    number_of_clusters = 5
    genes_for_clustering = ['Ppy', 'Ins1', 'Gcg', 'Ngp', 'S100a8', 'Cytl1', 'Pyy', 'Guca2a', 'Apoa1', 'Tff3', 'S100a9', 'Cldn5', 'Mgp', 'Hba-a1', 'Apoa2', 'Camp', 'Ctrb1', 'Mfap4', 'Iapp', 'C1qc', 'Cldn4', 'Prss2', 'Calm4', 'Nkg7', 'Try4', 'Lcn2', 'Chad', 'Penk', 'Lyz2', 'Krt14', 'Cela2a', 'H2-Aa', 'Aldoc', 'Cel', 'Ly6d', 'Sprr1a', 'Retnlb', 'C1qa', 'Ctrl', 'Phgr1', 'C1qb', 'Clps', 'Gpx3', 'Krt15', 'Sycn', 'Alb', 'H2-Ab1', 'Mir568', 'Csf3', '1810065E05Rik', 'Zg16', 'Cxcl10', 'Serpina3k', 'Mir744', 'Krt17', 'Cldn3', 'Ang4', 'Mup3', 'Gpr37l1', 'Gpihbp1', 'Cxcl1', 'Olig1', 'Cd74', 'Tacstd2', 'Aqp8', 'Agr2', 'Mup20', 'Rnase1', 'Ccl5', 'Fabp2', 'Myod1', 'Ccl2', 'Klk1', 'Krt6a', 'Sit1', 'Sval1', 'Ins2', 'Krt13', 'Ptprcap', 'Try5', 'Ptgds', 'Ccl7', 'Acta2', 'Nnat', 'Krt10', 'Igfbp6', 'Ly6c1', 'Krt4', 'Lgals4', 'Pnlip', 'Mir3061', 'Pi16', 'Gzma', 'Ahsg', 'Phxr4', 'Lgals7', 'Chga', 'Tmem119', 'Mt3', 'Tspan1', 'Hist1h2ae', 'Apoc3', 'Mir1934', 'Ctla2a', 'H2-Eb1', 'Krt5', 'Vtn', 'Mir3066', 'Aldob', 'Scg2', 'Cd8a', '2210010C04Rik', 'Calml3', 'Ctsw', 'Krt18', 'Cd68', 'Fabp1', 'Fgg', 'Il2rg', 'Snord73a', 'Cd7', 'Ly6c2', 'Sfn', 'Cldn7', 'Fahd1', 'Car1', 'H2-Q10', 'Hp', 'Cpa1', 'Reg1', 'Gzmb', 'Hilpda', 'Krtdap', 'Mir370', 'Cd14', 'Mt4', 'Fgb', 'Gimap9', 'Dcn', 'Mir24-2', 'Sox18', 'Mir326', 'Ctxn1', 'Mmp3', 'Spp1', 'Pnliprp1', 'Cd79b', 'Cela3b', 'Cyr61', 'AI467606', 'Ntsr2', 'Myoc', 'Gnmt', 'Apoe', 'Retnlg', 'Gpx2', 'Snord104', 'Clec4g', 'Ltf', 'Rnf186', 'Gpr182', 'Hist1h2ab', 'Ly6a', 'Resp18', '2200002D01Rik', 'Krt8', 'Bex2', 'Mir1931', 'Clec14a', 'Gpr84', 'Ccl4', 'Krt16', 'Apoc1', 'Lrg1', 'Snord47', 'Lum', 'S100a4', 'Mir221', 'Snord52', 'Mir142', 'Hspa1a', 'Ttr', 'Hes5', 'Tyrobp', 'Cd79a', 'Id1', 'S1pr1', 'Prss30', 'Hpx', 'Cela1', 'Ccl6', 'Apof', 'Serping1', 'Il1b', 'Krt20', 'Hspb1', 'Mxra8', 'Miox', 'H2-DMa', 'Hist2h2ac'] 
    doi = '00000'
    indexes = list(range(100)) 
    """
    
    #real code
    number_of_clusters = order['method']['parameters']['number_of_clusters']
    genes_for_clustering = order['method']['parameters']['genes_for_clustering']
    doi = set(order['basket']['doi']).pop()   # currently no cross datasets requests, use set() for simplicity
    indexes = order['basket']['indexes']

    dataset_composition, error_message, integrated_dim_red, candidate_TPM, candidate_genes = \
         merge_expression_matrix(order)
    if dataset_composition == 'single_dataset' or dataset_composition == 'cross_datasets':
        error = False
        error_message = ''
        cell_types = calculate_unsupervised_clustering(candidate_TPM, candidate_genes, \
            number_of_clusters, genes_for_clustering, method='louvain')
    elif dataset_composition == 'error':
        pass # if error, won't reached here.

    dataToClient = {

        "matrix":{
            "labels":{ 
                "cellTypes": cell_types
            }
        },

        "message": {
            "error": error,
            "error_message": error_message,
        },
         
    }

    return dataToClient

#========================================
#  general functions
#========================================

def search_gene_expression(order, current_user):

    """
    dataToServer= {

        "method": {
            "type": "search_gene_expression",
            "parameters": {
                "genes":["A", "B", "C"],
                "use_sum": 0
            }      
        },

        "basket":{
            "doi": ['00000', '00000', '00000', '00000'],
            "indexes": [0,1,2,3]
        },
    
    }
    """

    """
    # test code
    user_input_genes = ['Ppy', 'Ins1', 'Gcg']
    user_input_genes = ['Ppy'] 
    doi = '00000'
    indexes = list(range(100))
    use_sum = 0
    """
    
    # real code
    user_input_genes = order['method']['parameters']['genes']
    gene = user_input_genes[0]
    # GSE_number = set(order['basket']['doi']).pop()  # currently no cross datasets requests, use set() for simplicity
    indexes = order['basket']['indexes']
    use_sum = order['method']['parameters']['use_sum']

    #candidate_TPM, candidate_genes = query_cells_of_given_indexes(doi, indexes)
    #expression_vectors = query_gene_expression_value(candidate_TPM, candidate_genes, user_input_genes, use_sum)

    doi = order['basket']['doi']
    datasets = list(set(doi))
    if len(datasets) == 1:
        GSE_number = datasets[0]
        not_found, expression_vector = query_gene_expression_individual_dataset(GSE_number=GSE_number, gene=gene)
        expression_vectors = expression_vector
        error = False
        error_message = ''
        if not_found:
            error = True
            error_message = 'The gene is not found within this dataset'

    else:
        #if current_user.merged_matrix: # save for future use, with logged in user
        #    print('used matrix stored in session')
        
        dataset_composition, error_message, integrated_dim_red, candidate_TPM, candidate_genes = \
            merge_expression_matrix(order)
    
        if dataset_composition == 'cross_datasets':
            error = False
            error_message = ''
            try:
                expression_vector = np.log2(candidate_TPM[:, candidate_genes.tolist().index(gene)].ravel() + 1).tolist()
                expression_vectors = expression_vector
            except ValueError as e:
                expression_vector = []
                expression_vectors = expression_vector
                error = True
                error_message = 'Gene not found, or not identified as common gene between different datasets.'

        elif dataset_composition == 'error':
            pass # if error happens, users won't be able to trigger the gene search function.
        
    dataToClient = {
    
        "matrix":{
            "data": {
                "expression_vectors": expression_vectors
            },
        },

         "message": {
            "error": error,
            "error_message": error_message,
        },

    }

    return dataToClient


#####################################################################################################################################


#  request handling function


#####################################################################################################################################

@csrf_exempt
def index_all_datasets(request):
    """
    response = {
        'index_of_all_datasets': [
            {
                'title': 'RNA Sequencing of Single Human Islet Cells Reveals Type 2 Diabetes Genes',
                'abstract': '...'
                'GSE_number': 'GSE81608'
            },
            {
            },
        ]
    }
    """
    print(request.user)
    response = query_index_of_all_datasets()

    return JsonResponse(response)

@csrf_exempt
def explore_dataset_regex(request, GSE_number):
    """
    Parameters:
        GSE_number: int()
            e.g.: 00000
    """

    # currently GSE81608 works

    response = {'GSE': GSE_number}
    print(response)

    response = query_individual_dataset('GSE' + str(GSE_number))

    return JsonResponse(response)

@csrf_exempt
def explore_private_dataset_regex(request, GSE_number):
    """
    Parameters:
        GSE_number: int()
            e.g.: 00000
    """

    response = {'PRIVATE': GSE_number}
    print(response)

    response = query_individual_dataset('PRIVATE' + str(GSE_number))

    return JsonResponse(response)

@csrf_exempt
def explore_dataset_subtypes(request, GSE_number, cell_type):
    """
    Parameters:
        GSE_number: int()
            e.g.: 00000
        cell_type: str()
            e.g.: normal
    
    e.g.: GSE115978_normal
    """

    response = query_individual_dataset('GSE' + str(GSE_number) + '_' + cell_type)

    return JsonResponse(response)

@csrf_exempt
def query_gene_expression_all_cells(request):
    """
    request = {"GSE_number": "GSE11111", "gene": "Egfr"}
    """
    GSE_number = "GSE11111"
    gene = "Ins1"
    if request.method == 'POST':
        order = json.loads(request.body)
        GSE_number = order["GSE_number"]
        gene = order["gene"]
    
    tsne_1, tsne_2, _ = get_tsne(GSE_number)

    my_db_loader = db_loader(target_db = GSE_number).load_data_from_db()
    TPM, genes, cell_types = my_db_loader.matrix, my_db_loader.genes, my_db_loader.cell_types
    
    not_found = 'FALSE'
    try:
        expression_vector = np.log2(TPM[:, genes.tolist().index(gene)] + 1).tolist()
    except:
        expression_vector = []
        not_found = 'TRUE'

    response = {
        "tsne_1": tsne_1,
        "tsne_2": tsne_2,
        "expression_vector": expression_vector,
        "cell_types": cell_types.tolist(),
        "not_found": not_found,
    }

    return JsonResponse(response)

@csrf_exempt
def gene_hint(request):
    """
    Parameters:
    -----------
        request:
            request = {"GSE_number": '', "gene": ''}
    """

    GSE_number = "GSE109774_Bladder_10x"
    gene = "ol"
    if request.method == 'POST':
        order = json.loads(request.body)
        GSE_number = order["GSE_number"]
        gene = order["gene"]
    
    genes = get_genes(GSE_number)
    candidate_genes = [i for i in genes if re.search(gene, i, re.IGNORECASE)]

    if len(candidate_genes) > 10:
        candidate_genes = candidate_genes[:10]

    error = False 
    error_message = ''
    if len(candidate_genes) == 0:
        error = True
        error_message = 'no match found'

    response = {
        "candidate_genes": candidate_genes,
        "message": {
            "error": error,
            "error_message": error_message 
        }
    }

    return JsonResponse(response)

@csrf_exempt
def gene_hint_all_datasets(request):
    """
    Parameters:
    -----------
        request:
            request = {"gene": ''}
    """

    gene = 'Egfr'
    if request.method == 'POST':
        order = json.loads(request.body)
        gene = order["gene"]
    
    genes = get_all_genes_of_species(taxonomy_id='10090', type='symbol') \
                + get_all_genes_of_species(taxonomy_id='9606', type='symbol') \
                + get_all_genes_of_species(taxonomy_id='7955', type='symbol')
    candidate_genes = [i for i in genes if re.search(gene, i, re.IGNORECASE)]

    error = False 
    error_message = ''
    if len(candidate_genes) == 0:
        error = True
        error_message = 'no match found'

    response = {
        "candidate_genes": candidate_genes,
        "message": {
            "error": error,
            "error_message": error_message 
        }
    }

    return JsonResponse(response)

@csrf_exempt
def query_gene_expression_all_cells_new(request):
    """
    request = {"GSE_number": "GSE11111", "gene": "Egfr"}
    """
    GSE_number = "GSE0_ALIGNED_Mus_musculus_Bladder"
    gene = "Dcn"
    
    if request.method == 'POST':
        order = json.loads(request.body)
        GSE_number = order["GSE_number"]
        gene = order["gene"]
    
    tsne_1, tsne_2, cell_types = get_tsne(GSE_number)
    not_found, expression_vector = query_gene_expression_individual_dataset(GSE_number=GSE_number, gene=gene)

    error = False 
    error_message = ''
    if not_found:
        error = True
        error_message = 'The gene is not found within this dataset'

    response = {
        "tsne_1": tsne_1,
        "tsne_2": tsne_2,
        "cell_types": cell_types,
        "expression_vector": expression_vector,

        "message": {
            "error": error,
            "error_message": error_message 
        }

    }
    return JsonResponse(response)

@csrf_exempt
def submit_order(request):

    """
    SubmitOrderAjax   all detailed ajax should inherit this form

        $.ajax({

        url: 'http://118.190.148.166:8910/market/usr/cart_ajax',

        dataToServer= {

            // type of analysis
            "method": {
                "type": "",
                "parameters": {
                }      
            },
        
            "basket":{
                "basketID": [],
                "doi": ['00000', '00000', '00000', '00000'],
                "indexes": [0,1,2,3]
            },
        
        },
        
        dataToClient = {
        
            // versatile, defined explicitly below 
        
        }

    })

    ----------------------
    
    return a JsonResponse to a submit order request depending on the type of order
    """

    response = {}
    current_user = request.user
    if request.method == 'POST':
        order = json.loads(request.body)
        #print(order)
        order_type = order['method']['type']
        print(order_type)

        # return different responses according to the type of request

        # functions for unsupervised_analysis
        if order_type == 'plot_tsne':
            response = plot_tsne(order, current_user)

        elif order_type == 'unsupervised_clustering':
            response = unsupervised_clustering(order, current_user)

        # general functions
        elif order_type == 'search_gene_expression':
            response = search_gene_expression(order, current_user)

        # search_gene_across_all_datasets
        elif order_type == 'search_gene_across_all_datasets':
            response = search_gene_across_all_datasets(order)

    return JsonResponse(response)

@csrf_exempt
def scibet_choose_reference(request):
    """
    request = {
        "reference_GSE": str()
    },
    
    response = {

        "annotation":{
            {"x": array(), "y": array(), "cell_types": array()}
        },

        "projection_tsne":{
            {"tsne_1": array(), "tsne_2": array(), "cell_types": array()}
        },

    },
    """
    response = {}

    if request.method == 'POST':
        request = json.loads(request.body)
        GSE_number = request['reference_GSE']
        
        # return reference tsne
        tsne_1, tsne_2, cell_types, tsne_1_center, tsne_2_center, cell_types_center = get_tsne_center(GSE_number)
        
        response={

            "tsne_cor": {"tsne_1": tsne_1, "tsne_2": tsne_2, "cell_types": cell_types},

            "annotation":{"tsne_1_center": tsne_1_center, "tsne_2_center": tsne_2_center, "cell_types_center": cell_types_center},
        }

    return JsonResponse(response)

@csrf_exempt
def scibet_general(request):
    """
    request = {
        "reference_GSE": str()
    },
    
    """
    response = {}

    # identifying reference and test
    GSE_number = 'GSE11111'
    received_file = '/home/biodb/data/Abio/database/datasets_for_upload/mouse_test.csv'
    if request.method == 'POST':
        usr_request = request.POST
        print(usr_request)
        print(usr_request['reference_GSE'])
        GSE_number = usr_request['reference_GSE'] # use info from post
        received_file = request.FILES['file'] # use info from post

    # calculate_reference
    scibet_db_loader = db_loader(target_db=GSE_number)
    scibet_db_loader.load_data_from_db()
    my_classifier = ScibetClassifier()
    print(scibet_db_loader.cell_types)
    #my_classifier.calculate_core(scibet_db_loader.matrix, scibet_db_loader.genes, scibet_db_loader.cell_types)
    my_classifier.calculate_core(expr_matrix=scibet_db_loader.matrix, 
                                genes=scibet_db_loader.genes, 
                                cell_types=scibet_db_loader.cell_types,
                                select_genes=True, number_of_genes=1000, 
                                log_transformed=False)

    # processing test file    
    test_df = pd.read_csv(received_file, index_col=0)
    test_TPM = normalize_to_TPM(test_df.to_numpy())
    test_genes = test_df.columns.to_numpy()
    # predictions = my_classifier.predict_cell_type(test_TPM, test_genes).tolist()
    predictions = my_classifier.predict_cell_type(expr_matrix=test_TPM, 
                                genes=test_genes, log_transformed=False)
    predictions.tolist()

    # generate tsne cor
    reference_tsne_1, reference_tsne_2, reference_cell_types = get_tsne(GSE_number)
    projection_tsne_1, projection_tsne_2 = my_classifier.generate_projection(
        np.array(reference_tsne_1), np.array(reference_tsne_2), np.array(reference_cell_types)
    )
    
    # generate major cell_types
    from collections import Counter
    counts = Counter(predictions)
    major_cell_types = []
    for cell_type in counts:
        if counts[cell_type] / len(predictions) > 0.005:
            major_cell_types.append(cell_type)
    
    # generate marker genes for major cell_types
    markers = query_individual_dataset(GSE_number)['markers']
    reference_markers = [marker for marker in markers if marker['cell_type'] in major_cell_types]

    print('here ok')
    

    # generate heatmap:
    marker_genes_set = []
    marker_cell_types = []
    for marker in reference_markers:
        marker_genes_set += marker['marker_genes']
        marker_cell_types.append(marker['cell_type'])
    marker_genes_set = list(set(marker_genes_set))

    df_heatmap_selected_genes = my_classifier.df_TPM_averaged.reindex(columns=marker_genes_set, fill_value=0)
    df_heatmap = np.log1p(df_heatmap_selected_genes.loc[marker_cell_types,:])
    import seaborn as sns
    my_cluster_map = sns.clustermap(df_heatmap)
    row_index = my_cluster_map.dendrogram_row.reordered_ind
    column_index = my_cluster_map.dendrogram_col.reordered_ind
    df_clustermap = df_heatmap.iloc[row_index, column_index]
    heatmap_matrix = df_clustermap.to_numpy()

    # generate paga
    import scanpy as sc
    import anndata

    df_paga_prepare = pd.DataFrame(scibet_db_loader.matrix, scibet_db_loader.cell_types, scibet_db_loader.genes)
    if my_classifier.reference_cell_types.shape[0] > 10:
        df_paga = df_paga_prepare.loc[major_cell_types,:]
    else:
        df_paga = df_paga_prepare
    
    adata = anndata.AnnData(df_paga.to_numpy(), obs={"cell_types": df_paga.index}, var={"genes": df_paga.columns})
    sc.pp.recipe_zheng17(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
    sc.tl.diffmap(adata)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
    
    my_annotations = []
    cell_types_list = adata.obs['cell_types'].tolist()
    for i in range(0, len(cell_types_list)):
        my_annotations.append(cell_types_list[i])
    adata.obs["my_annotations"] = my_annotations
    sc.tl.paga(adata, groups='my_annotations')
    paga_info = adata.uns['paga']
    paga_nodes = adata.obs["my_annotations"].cat.categories.tolist()
    
    paga_connectivities = paga_info["connectivities"].toarray()
    my_connectivities = paga_connectivities
    unique = []
    for i in range(0, my_connectivities.shape[0]):
        for j in range(0, my_connectivities.shape[1]):
            if my_connectivities[i][j] != 0:
                if my_connectivities[i][j] not in unique:
                    unique.append(my_connectivities[i][j])
                else:
                    my_connectivities[i][j] = 0
    
    nodes = {
        "cell_type": [node for node in paga_nodes],
        "numver": [counts[node] for node in paga_nodes]
    }            
    my_connectivities = np.around(my_connectivities, 3).tolist()

    response = {

            'projection_plot':{
                'tsne_1': projection_tsne_1,
                'tsne_2': projection_tsne_2,
                'cell_types': predictions,
            },

            "major_cell_types": reference_markers,

            "heatmap": {
                "heatmap_cell_types": marker_cell_types, # array()
                "heatmap_genes": marker_genes_set, # array()
                "heatmap_matrix": heatmap_matrix.tolist() # nested array()
            },

            "paga": {
                "nodes" : nodes,
                "connectivities": my_connectivities,
            },

        }

    return JsonResponse(response)

@csrf_exempt
def scibet2(request):
    """
    response = {

            'tsne_plot':{
                'tsne_1': array(),
                'tsne_2': array(),
                'cell_types': array()
            },

            'projection_plot':{
                'tsne_1': array(),
                'tsne_2': array(),
                'cell_types': array(),
            }

        }    
    """

    if request.method == 'POST':

        received_FILE = request.FILES['file']
        received_file = pd.read_csv(received_FILE, index_col=0)
        print(received_FILE)
        print(type(received_FILE))
        
        received_file_name = received_FILE._get_name()
        print(received_file_name)

        # check file suffix
        if received_file_name.endswith('.csv'):

            train_core = pd.read_csv('./market/scibet/neuron_train_CORE.csv', index_col=0)
            test_x = received_file.iloc[:,1:]
            test_x=test_x.loc[:,train_core.columns].fillna(0)
            test_x1=np.log1p(test_x)
            pred_mtx=pd.DataFrame(np.dot(train_core,test_x1.T))
            pred_y=train_core.index[pred_mtx.idxmax(0)]
            predictions = pred_y.tolist()

            # assign tsne
            df_tsne = pd.read_csv('./market/scibet/brain_meta_sampled.csv')
            reference_tsne_1 = df_tsne['x'].tolist()
            reference_tsne_2 = df_tsne['y'].tolist()
            reference_cell_types = df_tsne['description'].tolist()

            projection_tsne_1 = []
            projection_tsne_2 = []

            for prediction in predictions:
                
                df_tsne_cell_type = df_tsne[df_tsne['description'] == prediction]
                tsne_1_list = df_tsne_cell_type['x'].tolist()
                tsne_2_list = df_tsne_cell_type['y'].tolist()
                try:
                    random_index = np.random.randint(0, len(tsne_1_list))
                    pseudo_tsne_1 = tsne_1_list[random_index]
                    pseudo_tsne_1 += float(np.random.rand(1) * 0.02 - 0.01)
                    projection_tsne_1.append(pseudo_tsne_1)
                    
                    random_index = np.random.randint(0, len(tsne_1_list))
                    pseudo_tsne_2 = tsne_2_list[random_index]
                    pseudo_tsne_2 += float(np.random.rand(1) * 0.02 - 0.01)
                    projection_tsne_2.append(pseudo_tsne_2)     
                except:
                    pass

    response = {

            #'tsne_plot':{
            #    'tsne_1': reference_tsne_1,
            #    'tsne_2': reference_tsne_2,
            #    'cell_types':reference_cell_types,
            #},

            'projection_plot':{
                'tsne_1': projection_tsne_1,
                'tsne_2': projection_tsne_2,
                'cell_types': predictions
            }

        }

    return JsonResponse(response)

    pass

@csrf_exempt
def scibet(request):

    response = {}

    # ==========================================================
    df_tsne = pd.read_csv('./market/scibet/tsne.csv')
    reference_tsne_1 = df_tsne['x'].to_numpy()
    reference_tsne_2 = df_tsne['y'].to_numpy()

    ident = df_tsne['ident']
    reference_cell_types = []
    core_cell_types = ['microglia','mesenchymal', 'B', 'endothelial', 'fibroblasts', 'T',
            'keratinocyte', 'epidermal', 'haematopoietic progenitor',
            'oligodendrocyte', 'epidermal', 'microglia', 'epithelial',
            'satellite cell', 'unknown', 'granulocytes', 'goblet', 'epithelial',
            'B', 'myeloid', 'epidermal', 'myeloid', 'epithelial',
            'endothelial', 'mesenchymal', 'endothelial', 'bladder cell',
            'urothelial', 'myeloid', 'T', 'non-beta islet', 'unknown', 'astrocyte',
            'beta islet', 'epithelial', 'hepatocytes', 'mesenchymal', 'NK',
            'neuron', 'mesenchymal', 'B', 'acinar',
            'oligodendrocyte precursors', 'granulocytes', 'epithelial',
            'ductal', 'endothelial', 'myeloid', 'endothelial',
            'cardiac muscle', 'T', 'epithelial', 'hepatocytes', 'unknown']

    for i in ident:
        reference_cell_types.append(core_cell_types[i])

    reference_cell_types = np.array(reference_cell_types)
    # ==========================================================
    cell_types_tsne = {}

    for i in range(0, len(core_cell_types)):
        _df = df_tsne[df_tsne['ident'] == i]
        tsne_1 = _df['tSNE_1'].tolist()
        tsne_2 = _df['tSNE_2'].tolist()
        cell_types_tsne[core_cell_types[i]] = {'tsne_1': tsne_1, 'tsne_2': tsne_2}
    # ==========================================================
    muris_train_core = pd.read_csv('./market/scibet/muris_train_core.csv')
    print('ok')
    # ==========================================================
    muris_test = pd.read_csv('./market/scibet/mouse_test.csv', index_col=0)
    received_file = muris_test

    # ==========================================================
    if request.method == 'POST':
        received_file = request.FILES['file']

        received_file_name = received_file._get_name()
        print(received_file_name)

        # check file suffix
        if received_file_name.endswith('.csv'):
            received_file = pd.read_csv(request.FILES['file'])  # the csv received from usr upload 
        
        else:
            pass

    # ==========================================================    
    reference_core_cell_types=np.array(core_cell_types)
    my_classifier = scibet_classifier()
    my_classifier.load_reference('./market/scibet/muris_train_core.csv', reference_core_cell_types, cell_types_tsne, reference_cell_types)

    test_df = received_file.loc[:, my_classifier.reference_genes.tolist()].fillna(0)
    test_matrix = test_df.to_numpy()
    my_classifier.predict_prob(test_matrix)

    # ========================================================== 
    projection_tsne_1 = []
    projection_tsne_2 = []

    for i in my_classifier.prediction_vector:
        group = core_cell_types.index(i)
        _df = df_tsne[df_tsne['ident'] == group]
        all_tsne_1 = _df['tSNE_1'].tolist()
        all_tsne_2 = _df['tSNE_2'].tolist()
        
        pseudo_tsne_1 = all_tsne_1[np.random.randint(0, len(all_tsne_1))]
        pseudo_tsne_1 += float(np.random.rand(1) * 0.2 - 0.1)
        projection_tsne_1.append(pseudo_tsne_1)

        pseudo_tsne_2 = all_tsne_2[np.random.randint(0, len(all_tsne_2))]
        pseudo_tsne_2 += float(np.random.rand(1) * 0.2 - 0.1)
        projection_tsne_2.append(pseudo_tsne_2)

    predicted_cell_types = my_classifier.prediction_vector.tolist()
    # =========================================================
    # =========================================================

    #raw_counts = muris_test.to_numpy()
    #genes = muris_test.columns.to_numpy()
    #TPM = normalize_to_TPM(raw_counts)
    #informative_genes, informative_t_scores, informative_TPM = select_informative_genes(TPM, genes, 500, labels = [])
    #unsupervised_tsne_1, unsupervised_tsne_2 = tsne(np.log2(TPM+1))
    
    response = {

            #'tsne_plot':{
            #    'tsne_1': unsupervised_tsne_1,
            #    'tsne_2': unsupervised_tsne_2,
            #},

            'projection_plot':{
                'tsne_1': projection_tsne_1,
                'tsne_2': projection_tsne_2,
                'cell_types': predicted_cell_types,
            }

        }

    return JsonResponse(response)

@csrf_exempt
def explore_dataset(request):
    """
    returns the basic info to explore a dataset

    response = {

        "info":{
            "title": '',
            "doi": '',
            "abstract": '',
        },

        "pre_calculated_tsne":{
            "tsne_1":[],
            "tsne_2":[],
            "cell_type":[]
        }

    }

    """
    
    response = query_individual_dataset('00000')

    return JsonResponse(response)

@csrf_exempt
def explore_usr_upload_dataset(request):
    """
    returns the basic info to explore a dataset

    response = {

        "pre_calculated_tsne":{
            "tsne_1":[],
            "tsne_2":[],
            "cell_type":[]
        },

        "markers": [
            {
                'cell_type': str(),
                'marker_genes': array(),
                "marker_t_scores": array()
            },
            {

            },
        ]

    }

    """

    response = {}
    
    if request.method == 'POST':

        # receive the file and load it into db
        # print('\n', request.FILES)
        received_file = request.FILES['file']
        print(received_file)
        print(type(received_file))
        
        received_file_name = received_file._get_name()
        print(received_file_name)

        # check file suffix
        if received_file_name.endswith('.csv'):
            df = pd.read_csv(received_file, index_col=0)
            matrix = df[df.columns[1:]].to_numpy()
            TPM = matrix / np.sum(matrix, axis=1).reshape(-1, 1) * 1e6
            cell_types = df['Unnamed: 0.1'].tolist()
            
        else:
            pass
        
        #=========================================
        target_db = 'test_upload'
        genes = df.columns[1:].to_numpy()
        usr_upload_db_loader = db_loader(genes, np.array(cell_types), TPM, target_db)
        usr_upload_db_loader.load_matrix()

        from sklearn.decomposition import PCA
        result = PCA(n_components=2).fit_transform(np.log2(TPM + 1))

        tsne_1, tsne_2 = result[:,0], result[:,1]
        #=========================================

        genes = df.columns[1:].to_numpy()
        markers = []

        for cell_type in set(cell_types):
            df_cell_type = df[df['Unnamed: 0.1'] == cell_type]
            matrix_cell_type = df_cell_type[df_cell_type.columns[1:]].to_numpy()
            TPM_cell_type = matrix_cell_type / np.sum(matrix_cell_type, axis=1).reshape(-1, 1) * 1e6
            informative_genes, informative_t_scores, informative_TPM, p_vals = \
                select_informative_genes(TPM_cell_type, genes, 50)
            marker = {
                'cell_type': cell_type,
                'marker_genes': informative_genes.tolist(),
                "marker_t_scores": informative_t_scores.tolist()
            }
            markers.append(marker)

        response = {

            "pre_calculated_tsne":{
                "tsne_1":tsne_1.tolist(),
                "tsne_2":tsne_2.tolist(),
                "cell_type":cell_types,
            },

            "markers": markers,

            'new_GSE': 'test_upload'
            
        }

        print(response)

    return JsonResponse(response)