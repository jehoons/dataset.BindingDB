import pandas as pd 
import os,sys,csv,gzip 
from os.path import join,dirname,exists 
from ipdb import set_trace 
from os import system 

import pprint
import collections
import operator

import pandas
import requests
import pickle

from downloader import download


scratch_dir = join(os.getcwd(), '/scratch')

datasets_url = 'http://192.168.0.97/share/StandigmDB/datasets'

url_dataset_GeneID = join(datasets_url, 
    'binding_db/GeneID.tsv.gz'
    )

url_dataset_BindingDBAll = join(datasets_url, 
    'binding_db/2015-10/BindingDB_All_2015m10.tsv'
    )

target_fields = [
    'BindingDB Target Chain  Sequence',
    'PDB ID(s) of Target Chain',
    'UniProt (SwissProt) Recommended Name of Target Chain',
    'UniProt (SwissProt) Entry Name of Target Chain',
    'UniProt (SwissProt) Primary ID of Target Chain',
    'UniProt (SwissProt) Secondary ID(s) of Target Chain',
    'UniProt (SwissProt) Alternative ID(s) of Target Chain',
    'UniProt (TrEMBL) Submitted Name of Target Chain',
    'UniProt (TrEMBL) Entry Name of Target Chain',
    'UniProt (TrEMBL) Primary ID of Target Chain',
    'UniProt (TrEMBL) Secondary ID(s) of Target Chain',
    'UniProt (TrEMBL) Alternative ID(s) of Target Chain',
    ]

chains_key = 'Number of Protein Chains in Target (>1 implies a multichain complex)'

chk_bindingdb_step1_pkl = join(scratch_dir, 
    'chk-bindingdb-step1.pkl'
    )

out_bindingdb_dataset = join(scratch_dir, 
    'processed-bindingdb-dataset.tsv'
    )


def read_bindingdb(path, verbose=False, max_rows=None):
    """
    Field documentation: https://www.bindingdb.org/bind/chemsearch/marvin/BindingDB-TSV-Format.pdf
    """
    reader = csv.reader(open(path,'rt'), delimiter='\t')
    header = next(reader)
    chains_index = header.index(chains_key)
    
    target0_index = chains_index + 1
    ligand_fields = header[:chains_index + 1]
    for j, row in enumerate(reader):
        if max_rows is not None and j == max_rows:
            break
        row = [x if x else None for x in row]
        ligand_values = row[:chains_index + 1]
        # Ensure line has sufficient ligand fields
        if len(row) < chains_index + 1:
            if verbose:
                print('Line', j + 2, 'is deficient')
            continue
        rowdict = collections.OrderedDict(zip(ligand_fields, ligand_values))
        for key in [chains_key]:
            if key not in rowdict:
                print(j+2)
                print(row)
                print(rowdict)
            rowdict[key] = int(rowdict[key])
        chains = list()
        assert rowdict[chains_key] == len(row[target0_index:]) / len(target_fields)
        for i in range(rowdict[chains_key]):
            i_0 = target0_index + i * len(target_fields)
            i_1 = target0_index + (i + 1) * len(target_fields)
            target_values = row[i_0:i_1]
            chain = collections.OrderedDict(zip(target_fields, target_values))
            chains.append(chain)
        rowdict['chains'] = chains
        yield rowdict


def step1(force=False):

    if exists(chk_bindingdb_step1_pkl) and force==False: 
        with open(chk_bindingdb_step1_pkl, 'rb') as f: 
            uniprot_to_entrez = pickle.load(f)

    else: 
        local_filename = download(remote_path=url_dataset_GeneID, 
            download_dir=scratch_dir, force=True
            )
        uniprot_df = pandas.read_table(local_filename, 
            compression='gzip'
            )
        uniprot_to_entrez = dict()
        
        for uniprot, entrez in zip(uniprot_df.uniprot, uniprot_df.GeneID):
            uniprot_to_entrez.setdefault(uniprot, set()).add(str(entrez))

        with open(chk_bindingdb_step1_pkl, 'wb') as f: 
            pickle.dump(uniprot_to_entrez, f)

    return uniprot_to_entrez


def step2():

    local_filename = download(remote_path=url_dataset_BindingDBAll,
        download_dir=scratch_dir)

    return local_filename    


def step3(mapdata=None, bindingdbpath=None):

    bindingdb_generator = read_bindingdb(bindingdbpath, verbose=True)
    bindings = list()
    
    # from tqdm import tqdm 
    
    for row in bindingdb_generator:
        if len(row['chains']) != 1:
            continue

        chain, = row['chains']
        uniprots = chain['UniProt (SwissProt) Primary ID of Target Chain']
        if not uniprots:
            continue
        
        uniprots = uniprots.split(',')

        template = dict()
        template['bindingdb_id'] = row['BindingDB MonomerID']
        template['reaction_id'] = row['BindingDB Reactant_set_id']
        template['source'] = row['Curation/DataSource']
        template['organism'] = row['Target Source Organism According to Curator or DataSource']
        template['pubmed'] = row['PMID']
        template['doi'] = row['Article DOI']

        affinities = {'Ki': row['Ki (nM)'], 'Kd': row['Kd (nM)'], 'IC50': row['IC50 (nM)']}
        for measure, affinity in affinities.items():
            if affinity is None:
                continue
            for uniprot in uniprots:
                entrez_set = mapdata.get(uniprot)
                if not entrez_set:
                    # uniprot_id not found in mapping
                    continue
                for entrez in entrez_set:
                    binding = template.copy()
                    binding['measure'] = measure
                    binding['affinity_nM'] = affinity
                    binding['uniprot'] = uniprot
                    binding['entrez_gene'] = entrez
                    bindings.append(binding)

    # Convert affinities to floats
    lt, gt, eq, err = 0, 0, 0, 0
    for binding in bindings:
        affinity = binding['affinity_nM']
        if affinity.startswith('<'):
            affinity = affinity.lstrip('<')
            affinity = float(affinity)
            if affinity >= 10.0:
                affinity -= 1.0
            lt += 1
        elif affinity.startswith('>'):
            affinity = affinity.lstrip('>')
            affinity = float(affinity)
            affinity += 1.0
            gt += 1
        else:
            try:
                affinity = float(affinity)
                eq += 1
            except ValueError:
                affinity = None
                err += 1

        binding['affinity_nM'] = affinity
    
    print('< {}\n> {}\n= {}\nerrors {}'.format(lt, gt, eq, err))

    fields = ['reaction_id', 'bindingdb_id', 'uniprot', 'entrez_gene',
              'measure', 'affinity_nM', 'source', 'organism', 'pubmed', 'doi']

    with open(out_bindingdb_dataset, 'wt') as write_file:
        writer = csv.DictWriter(
            write_file, 
            delimiter='\t', 
            fieldnames=fields
            )

        writer.writeheader()
        
        bindings.sort(
            key=operator.itemgetter(*fields)
            )
        
        writer.writerows(bindings)


    # extra information: 

    # Measurement types
    # measure_keys = ['Ki (nM)', 'IC50 (nM)', 'Kd (nM)', 'EC50 (nM)'] #, 'kon (M-1-s-1)', 'koff (s-1)']
    # bindingdb_generator = read_bindingdb(step2())
    # measures = list()
    # for i, row in enumerate(bindingdb_generator):
    #     if len(row['chains']) != 1:
    #         continue
    #     chain, = row['chains']
    #     uniprot = chain['UniProt (SwissProt) Primary ID of Target Chain']
    #     if not uniprot:
    #         continue
    #     measure_set = frozenset(key for key in measure_keys if row[key] is not None)
    #     measures.append(measure_set)

    # pprint.pprint(collections.Counter(measures))

    # # Number of chains (proteins in target)
    # bindingdb_generator = read_bindingdb(step2())
    
    # collections.Counter(int(row[chains_key]) 
    #     for row in bindingdb_generator)

    # # Targets that mapped to SwissProt
    # bindingdb_generator = read_bindingdb(step2())

    # collections.Counter(
    #     bool(row['chains'][0]['UniProt (SwissProt) Primary ID of Target Chain'])
    #     for row in bindingdb_generator if len(row['chains']) == 1
    # )

    # # Species
    # bindingdb_generator = read_bindingdb(step2())
    # collections.Counter(
    #     row['Target Source Organism According to Curator or DataSource']
    #     for row in bindingdb_generator if
    #     len(row['chains']) == 1 and 
    #     row['chains'][0]['UniProt (SwissProt) Primary ID of Target Chain']
    # ) 


def test_run():

    uniprot_to_entrez = step1()
    
    local_bindingdb_path = step2()

    if not exists(out_bindingdb_dataset): 
        step3(mapdata=uniprot_to_entrez, 
            bindingdbpath=local_bindingdb_path)

    assert True 






