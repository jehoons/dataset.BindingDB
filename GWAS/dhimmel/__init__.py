# import collections

from ipdb import set_trace
from os.path import exists,dirname,join,basename 
import pandas as pd, numpy as np  

module_dir = dirname(__file__)


def load_gene_assoc():
    filepath = join(module_dir, 'data', 
        'gene-associations.tsv'
        )
    
    return pd.read_csv(filepath, sep='\t')


def load_snp_assoc():
    filepath = join(module_dir, 'data', 
        'snp-associations.tsv'
        )

    return pd.read_csv(filepath, sep='\t')


def test_module():
    from hetio.io import hetiodb
    from hetio import anal

    dataset_DaG = load_gene_assoc()[['doid_name', 'symbol']]
    # data2 = load_snp_assoc()

    clone_links = dataset_DaG

    clone_links_set = set(['%s-%s'%(x[0].strip(),x[1].strip()) 
        for x in clone_links.values.tolist()])

    hetio_links = hetiodb.return_relation_names(rtype='ASSOCIATES_DaG')    
    hetio_links_set = set(['%s-%s'%(x[0].strip(),x[1].strip()) 
        for x in hetio_links.values.tolist()]) 

    res = anal.inclusion_analysis(hetio_data=hetio_links_set, \
        clone_data=clone_links_set)

    from pprint import pprint 

    pprint(res)

    set_trace()
