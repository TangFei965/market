######################################

#regarding db design:
#1. follow the Anndata style
#2. use pickle for compression
#3. store raw counts and then do the calculations
#4. gene t scores should be rounded 
#5. TPM should also be rounded

######################################

######################################

######################################

from pymongo import MongoClient

from .algorithms import *

class db_loader():
    
    def __init__(self, genes=None, cell_types=None, matrix=None, target_db=None):
        """
        Parameters:
        -----------
        genes: np.array()
            genes names
            n genes
            
        cell_types:np.array()
            cell_type
            m cells
        
        matrix: np.array()
            might be raw counts or TPM, as long as not logged.
            will transform into TPM before calculations
            m by n matrix
            
        target_db: str()
            the target db to insert into
        """
        self.genes = genes
        self.cell_types = cell_types
        self.matrix = matrix
        self.target_db = target_db

    def load_data_from_db(self):
        """
        load the basic info from the matrix and the annotations
        """

        client = MongoClient()
        db = client[self.target_db]
        
        # get genes
        gene_collection = db['gene']
        genes = np.array(gene_collection.find_one({"type": "symbol"})['content'])
        
        # get cell_types
        cell_type_collection = db['pre_calculated_tsne']
        cell_types = np.array(cell_type_collection.find_one({'field':'cell_type'})['content'])

        # get expression matrix
        matrix_collection = db['matrix']
        number_of_cells = cell_types.shape[0]
        indexes = list(range(0, number_of_cells))
        rows = []
        for index in indexes:
            #print('doing %s of %s' %(index, number_of_cells))
            row = matrix_collection.find_one({"index": index})['matrix']['data']
            rows.append(row)         
        matrix = np.vstack(rows)
        
        self.genes = genes
        self.cell_types = cell_types
        self.matrix = matrix

        return self

    def load_info(self, title='empty', abstract='empty', doi='empty', doi_link='empty', figure_url='empty',):
        
        # rewrite using a for loop
        """
        schema:
        
            {"field": "title", "content": ""},
            {"field": "abstract", "content": ""},
            {"field": "doi", "content": ""},
            {"field": "doi_link", "content": ""},
            {"field": "figure_url", "content": ""},

            # regarding experimental design and methodology
            {"field": "platform"}

            {"field": "species", "content": ""},
            {"field": "tissue", "content": ""}, # tissue name should follow http://gepia.cancer-pku.cn/detail.php?gene=
            {"field": "disease", "content": ""},
        """
        
        client = MongoClient()
        db = client[self.target_db]  # self replace
        collection = db['info']
        
        def _load_title(_title):
            try:
                existing_title = db['info'].find({"field": "title"}).next()['content']
                print('\n', 'title: \n\t %s \n has exist, please consider whether replace it' % existing_title)
                print("if you see the content is 'empty', it is to keep database schema consistency, safely add a real content")
            except StopIteration:
                document = {}
                document['field'] = 'title'
                document['content'] = _title
                collection.insert_one(document)
                print('title: \n %s \n is loaded successfully' % _title)
                
            print('-----------------------------------')
        
        def _load_abstract(_abstract):
            try:
                existing_abstract = db['info'].find({"field": "abstract"}).next()['content']
                print('\n', 'abstract: \n\t %s \n has exist, please consider whether replace it' % existing_abstract)
                print("\nif you see the content is 'empty', it is to keep database schema consistency, safely add a real content")
            except StopIteration:
                document = {}
                document['field'] = 'abstract'
                document['content'] = _abstract
                collection.insert_one(document)
                print('abstract: \n %s \n is loaded successfully' % _abstract)
            
            print('-----------------------------------')
        
        def _load_doi(_doi):
            try:
                existing_doi = db['info'].find({"field": "doi"}).next()['content']
                print('\n', 'doi: \n\t %s \n has exist, please consider whether replace it' % existing_doi)
                print("\nif you see the content is 'empty', it is to keep database schema consistency, safely add a real content")
            except StopIteration:
                document = {}
                document['field'] = 'doi'
                document['content'] = _doi
                collection.insert_one(document)
                print('doi: \n %s \n is loaded successfully' % _doi)
                
            print('-----------------------------------')
        
        def _load_doi_link(_doi_link):
            try:
                existing_doi_link = db['info'].find({"field": "doi_link"}).next()['content']
                print('\n', 'doi_link: \n\t %s \n has exist, please consider whether replace it' % existing_doi_link)
                print("\nif you see the content is 'empty', it is to keep database schema consistency, safely add a real content")
            except StopIteration:
                document = {}
                document['field'] = 'doi_link'
                document['content'] = _doi_link
                collection.insert_one(document)
                print('doi_link: \n %s \n is loaded successfully' % _doi_link)
                
            print('-----------------------------------')
        
        def _load_figure_url(_figure_url):
            try:
                existing_figure_url = db['info'].find({"field": "_figure_url"}).next()['content']
                print('\n', 'figure_url: \n\t %s \n has exist, please consider whether replace it' % existing_figure_url)
                print("\nif you see the content is 'empty', it is to keep database schema consistency, safely add a real content")
            except StopIteration:
                document = {}
                document['field'] = '_figure_url'
                document['content'] = _figure_url
                collection.insert_one(document)
                print('_figure_url: \n %s \n is loaded successfully' % _figure_url)
                
            print('-----------------------------------')
            
        _load_title(title)
        _load_abstract(abstract)
        _load_doi(doi)
        _load_doi_link(doi_link)
        _load_figure_url(figure_url)

    def load_matrix(self):
        """
        
        final documents in the matrix collection
        
        {"genes": ['a', 'b', 'c']},

            x n      
                { 
                    "_id" : ObjectId("542c2b97bac0595474108b48"),  system generated
                    
                    "index": 0,

                    "matrix": {

                        "labels":{
                            "cell_type": ...
                        }

                        "data":[] need to be changed into sparse matrix, but csr is not json serializable...
                    }
                
                }
        
        """
        
        # preparing documents
        documents = [{"genes": self.genes.tolist()}]

        for i in range(0, len(self.cell_types)):
            cell = {"index": i}
    
            matrix = {}
            matrix["labels"] = {"cell_type": self.cell_types.tolist()[i]}
    
            m = self.matrix[i,:].tolist()
            matrix["data"] = m
    
            cell["matrix"] = matrix
    
            documents.append(cell)
        
        # load documents into db
        client = MongoClient()
        db = client[self.target_db]
        db.drop_collection('matrix')
        collection = db['matrix']
        collection.insert_many(documents)

        return 'success!'

    def load_transposed_matrix(self):
        """
        
        gene by cell matrix to accelerate reading speed.
        
        final documents in the matrix collection
        
        {
            x n      
                { 
                    "_id" : ObjectId("542c2b97bac0595474108b48"),  system generated
                    
                    "index": 0,
                    
                    "gene": 'abc',

                    "data":[],
                    }
                
                }
        
        """
        
        # get gene by cell matrix
        transposed_matrix = self.matrix.T
        
        # preparing documents
        documents = []
        for i in range(0, self.genes.shape[0]):
            current_gene = {"index": i}
            current_gene['gene'] = self.genes[i]
            current_gene['data'] = transposed_matrix[i, :].tolist()   
            documents.append(current_gene)
        
        # load documents into db
        client = MongoClient()
        db = client[self.target_db]
        collection = 'transposed_matrix'
        db.drop_collection(collection)
        collection = db[collection]
        collection.insert_many(documents)

        return 'success!'

    def load_pre_calculated_tsne(self):
        # note this method now has an intrinsic problem, that a tsne map should be fine tuned manually before loading into db.
        # re-design it in the following way:
            # load matrix from db, calculate tsne with parameters free
            # add a confirm load tsne helper function
        """
        final documents in the pre_calculated_tsne collection
        
        {
        "field": "tsne_1",
        "content": list()
        },
        
        {
        "field": "tsne_2",
        "content": list()
        },
        
        {
        "field": "cell_type",
        "content": list()
        }
        """
        
        # calculate tsne 
        TPM = normalize_to_TPM(self.matrix)
        informative_genes, informative_t_scores, informative_TPM, _ = \
            select_informative_genes(TPM, self.genes, 500, labels = [])
            
        tsne_1, tsne_2 = tsne(np.log2(informative_TPM + 1))
        
        # preparing documents
        documents = [
            {
                "field": "tsne_1",
                "content": tsne_1
            },
            {
                "field": "tsne_2",
                "content": tsne_2
            },
            {
                "field": "cell_type",
                "content": self.cell_types.tolist()
            }
        ]
        
        # load documents into db
        client = MongoClient()
        db = client[self.target_db]
        db.drop_collection('pre_calculated_tsne')
        collection = db['pre_calculated_tsne']
        collection.insert_many(documents)

        return tsne_1, tsne_2 
    
    def load_manually_calibrated_tsne(self, tsne_1, tsne_2):
        # consider combine with the one above into one function using manual calibration = False
        # input should be arrays
        """
        final documents in the pre_calculated_tsne collection
            
        {
        "field": "tsne_1",
        "content": list()
        },
            
        {
        "field": "tsne_2",
        "content": list()
        },
            
        {
        "field": "cell_type",
        "content": list()
        }
        
        Parameters:
            tsne_1: np.array()
            tsne_2: np.array()
        ----------------
        """
            
        client = MongoClient()
        db = client[self.target_db]
        collection = db['pre_calculated_tsne']
        
        if collection.count_documents({}) != 0:
            warning_message = 'this collection has already has data, please decide whether to rewrite it'
            return warning_message
        
        documents = [
            {
                "field": "tsne_1",
                "content": tsne_1.tolist()
            },
            {
                "field": "tsne_2",
                "content": tsne_2.tolist()
            },
            {
                "field": "cell_type",
                "content": self.cell_types.tolist()
            }
        ]
        
        collection.insert_many(documents)
        
        return 'success'

    def load_marker_genes(self, number_of_genes=50):
        
        client = MongoClient()
        db = client[self.target_db]
        db.drop_collection('marker_genes')
        collection = db['marker_genes']
        
        labels = self.cell_types.tolist()
        cell_types = set(self.cell_types.tolist())

        markers = []

        for cell_type in cell_types:
            modified_labels = []
            for i in range(0, len(labels)):
                if labels[i] == cell_type:
                    modified_labels.append(cell_type)
                else:
                    modified_labels.append('AVERAGE')
            
            marker_genes, marker_t_scores, informative_TPM, p_vals = \
                select_informative_genes(self.matrix, self.genes, number_of_genes, labels = modified_labels)
            
            markers.append([cell_type, marker_genes.tolist(), marker_t_scores.tolist(), p_vals.tolist()])

            document = {
                "cell_type": cell_type,
                "marker_genes": marker_genes.tolist(),
                "marker_t_scores": marker_t_scores.tolist(),
                "p_vals": p_vals.tolist(),
            }
            collection.insert_one(document)

        return markers

    def load_scibet_core(self):
        """
        DB_schema:
        -----------
            {"reference_genes": reference_genes:list()}:
                    number_of_documents = 1
                    
            {"cell_type": cell_type, "prob_vector": vector.tolist()}:
                    number_of_documents = number_of_cell_types
                    prob_vector: 
                        log2(X) transformed probability vector. 
                        X is the probability matrix normalized to sum = 1
        
        Parameters:
        -----------
            self
            basically, using self.matrix, self.genes, self.cell_types to calculate scibet core
        
        Returns:
        -----------
            loading_message: str()
        """
        
        # build train core
        my_classifier = ScibetClassifier()
        my_classifier.calculate_core(expr_matrix=self.matrix, genes=self.genes, 
                                    cell_types=self.cell_types, log_transformed=False,
                                    select_genes=True, number_of_genes=1000)
        reference_core = my_classifier.reference_core
        reference_genes = my_classifier.reference_genes
        reference_cell_types = my_classifier.reference_cell_types
        
        # preparing documents
        documents = [
            {"reference_genes": reference_genes.tolist()},
        ]
        for i in range(reference_cell_types.shape[0]):
            cell_type = reference_cell_types[i]
            vector = reference_core[i, :]
            documents.append({"cell_type": cell_type, "prob_vector": vector.tolist()})
        
        # write to db
        client = MongoClient()
        db = client[self.target_db]
        collection_name = 'scibet_core'
        db.drop_collection(collection_name)
        collection = db[collection_name]
        collection.insert_many(documents)
        
        message = 'finished loading scibet core, \n identified %d cell types' %(len(documents) - 1)
        print(message)

        return message
        
    
    
    