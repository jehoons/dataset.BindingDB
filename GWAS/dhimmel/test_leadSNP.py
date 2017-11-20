import collections

from ipdb import set_trace
from os.path import exists,dirname,join,basename 
import pandas as pd, numpy as np  
import pickle,json, requests
import os,math,re

from hetio.datasets import datasets_url
from hetio.preproc.downloader import download 

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn


test_dir = dirname(__file__)

download_dir = join(test_dir, 'download')

# Read EBI's GWAS Catalog with ontology annotations
# https://www.ebi.ac.uk/gwas/docs/fileheaders
url_gwas_catalog = join(datasets_url, 'gwas/dhimmel', 
    'gwas_catalog_v1.0.1-downloaded_2015-06-08.tsv.gz'
    )

url_symbols_human = join(datasets_url, 'gwas/dhimmel', 
    'entrez-symbols-human.json'
    )

# Read entrez synonym to genes mapping
url_synonyms_human = join(datasets_url, 'gwas/dhimmel',  
    'entrez-synonyms-human.json'
    )

# Read entrez symbol to gene mapping
url_genes_human = join(datasets_url, 'gwas/dhimmel', 
    'entrez-genes-human.tsv'
    )

# Read DO Slim propagated cross-references 
url_xrefs_prop_slim = join(datasets_url, 'gwas/dhimmel', 
    'disease-ontology-xrefs-prop-slim.tsv'
    )

# see_ebi_df = join(test_dir, 'data/see-ebi-df.html')

outfile_ebi = join(test_dir, 'data/snap/ebi-df.tsv')

outfile_lead_snp = join(test_dir, 'data/snap/do-slim-lead-SNPs.txt')


def _tad(df):
    # open(see_ebi_df,'w').write(ebi_df.head().to_html())
    from tempfile import NamedTemporaryFile
    myf = NamedTemporaryFile(delete=False)
    df.head().to_csv(myf.name); 
    os.system('tad -f {}'.format(myf.name))


def _preproc(): 
    dtype_dict = {
        'UPSTREAM_GENE_ID': np.character,
        'DOWNSTREAM_GENE_ID': np.character,
        'SNP_GENE_IDS': np.character
        }
    
    ebi_df = pd.read_csv(url_gwas_catalog, compression='gzip', 
        low_memory=False, dtype=dtype_dict, sep='\t') 

    # tad(ebi_df)

    #   Entrez Genee    
    localfile = download(
        remote_path=url_symbols_human, 
        download_dir=download_dir
        )
    symbol_to_id = json.load(open(localfile))

    localfile = download(remote_path=url_synonyms_human, 
        download_dir=download_dir
        )
    synonym_to_ids = json.load(open(localfile))

    # Read entrez genes
    localfile = download(
        remote_path=url_genes_human, download_dir=download_dir
        )
    
    entrez_df = pd.read_table(localfile)

    ## Create a set of coding genes
    coding_genes = set(entrez_df.GeneID[
        entrez_df.type_of_gene == 'protein-coding'
            ].astype(str))
  
    num_coding_genes = len(coding_genes)
    
    num_genes = len(entrez_df)

    ## Create a symbol dataframe
    symbol_df = entrez_df[['GeneID', 'Symbol']].rename(
        columns={'GeneID': 'gene', 'Symbol': 'symbol'}
        )
    
    symbol_df.gene = symbol_df.gene.astype(str)

    #
    #   Convert from EFO to Disease Ontology
    #

    ## Create a uri_df (for cross-references)
    rows = list()
    for uri in filter(pd.notnull, set(ebi_df['MAPPED_TRAIT_URI'])):
        head, tail = uri.rsplit('/', 1)
        resource, resource_id = tail.split('_', 1)
        rows.append([uri, resource, resource_id])
        
    uri_df = pd.DataFrame(rows, 
        columns=['MAPPED_TRAIT_URI', 'resource', 'resource_id'])

    localfile = download(
        remote_path=url_xrefs_prop_slim, 
        download_dir=download_dir
        )
    
    doxref_df = pd.read_table(localfile)

    # Inner join the GWAS catalog with the DO slim mapping 
    map_df = uri_df.merge(doxref_df)
    map_df = map_df[['MAPPED_TRAIT_URI', 'doid_code', 'doid_name']]
    ebi_df = ebi_df.merge(map_df)

    num_ebi_df = len(ebi_df)

    #
    #   Identify lead SNPs and proxy search using SNAP
    #

    # The SNAP proxy search must be run MANUALLY. Windows are computed in a 
    # companion notebook. Write all lead SNPs for SNAP input
    gwas_snps = set('rs{}'.format(x) 
        for x in ebi_df.SNP_ID_CURRENT if pd.notnull(x)
        )

    gwas_snps = sorted(gwas_snps)
    
    os.system('mkdir -p {}'.format(dirname(outfile_lead_snp)))
    with open(outfile_lead_snp, 'w') as write_file:
        write_file.write('\n'.join(gwas_snps))
 
    ebi_df.to_csv(outfile_ebi, index=False, sep='\t')

    return ebi_df, gwas_snps


def load():
    if exists(outfile_ebi) and exists(outfile_lead_snp): 
        ebi_df = pd.read_csv(outfile_ebi, sep='\t') 

        with open(outfile_lead_snp) as read_file:
            leads = {line.rstrip() for line in read_file}

    else:
        ebi_df, leads = _preproc()

    return {
        'GWAS-Catalog': ebi_df, 
        'LeadSNP': leads
        }


def test_load():

    dataset = load() 

    # set_trace()


