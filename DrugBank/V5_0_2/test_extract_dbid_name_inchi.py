import pandas as pd 
import os
import sys 
from os.path import join,dirname,exists
from ipdb import set_trace
from os import system
from hetio.downloader import download_v2
from hetio.datasets import datasets_url
import pytest
import json, gzip
from xmljson import badgerfish as bf, abdera, cobra, gdata, parker, yahoo
import xml.etree.ElementTree as ET
from collections import OrderedDict

import pybel 
from pybel import *

from hetio import scratch_dir
from hetio.preproc.DrugBank.V5_0_2 import loader as db_loader 

data_dir = dirname(__file__) + '/data'
os.system('mkdir -p {}'.format(data_dir))

chk_dbid_name_inchi = data_dir+'/chk-dbid-name-inchi.tsv'


def test_main():
    if exists(chk_dbid_name_inchi) : 
        return

    dataset = db_loader.load_table() 
    import numpy as np
    from tqdm import tqdm

    drugbank_inchi_list = [] 
    for idx in tqdm(dataset.index): 
        name = dataset.loc[idx, 'name']
        dbid = dataset.loc[idx, 'drugbank-id']
        dbid = eval(dbid)

        if type(dbid) ==list: 
            dbid = dbid[0]['content']
        elif type(dbid) == dict: 
            dbid = dbid['content']
        else: 
            assert False 

        cal = dataset.loc[idx, 'calculated-properties']

        record_template = {
            'name': name, 
            'DBID': dbid, 
            'InChI': None, 
            'InChIKey': None, 
            }

        if pd.isnull(cal):
            cal = '{}'

        caldict = eval(cal)
        if caldict == {}:
            drugbank_inchi_list.append(record_template)
            continue

        property = caldict['property']
        if type(property) != list:
            property = [property]

        df0 = pd.DataFrame(property).set_index('kind')
        
        if 'InChI' not in df0.index: 
            drugbank_inchi_list.append(record_template)
            continue 

        # set_trace()
        record_template['InChI'] = df0.loc['InChI', 'value']
        # set_trace()
        record_template['InChIKey'] = df0.loc['InChIKey', 'value']
        drugbank_inchi_list.append(record_template)

    df2 = pd.DataFrame(drugbank_inchi_list)
    df2.to_csv(chk_dbid_name_inchi, index=False, sep='\t')

    # set_trace()
