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

# license = ' CC BY-NC 4.0'

# base_dir = dirname(__file__)

data_dir = dirname(__file__) + '/data'; 
os.system('mkdir -p {}'.format(data_dir))

#
# input files 
#

url_drugbank_dataset = join(datasets_url, 'drug_bank/5-0-2',
    'drug-bank-full-5.0.2.xml.gz'
    )

url_drugbank_structure = join(datasets_url, 'drug_bank/5-0-2', 
    'structures.sdf'
    )

#
# output files 
#

# intermediate check files: 
chk_drugbank_1 = join(data_dir, 'chk-drugbank-1.json')
chk_drugbank_2 = join(data_dir, 'chk-drugbank-2.json')
chk_drugbank_3 = join(data_dir, 'chk-drugbank-3.tsv')
chk_drugbank_4 = join(data_dir, 'chk-drugbank-4(small).tsv')


def process_step1(force=False):
    if exists(chk_drugbank_1) and force==False:
        return 

    local_filename = download_v2(
        remote_path=url_drugbank_dataset, 
        download_dir=data_dir
        )

    with gzip.open(local_filename, 'rb') as f:
        xmldata = f.read()
    
    dictorized = yahoo.data(ET.fromstring(xmldata))

    with open(chk_drugbank_1, 'w') as fp:
        json.dump(dictorized, fp, indent=4)


def process_step2(force=False):
    if exists(chk_drugbank_3) and force==False:
        return 

    with open(chk_drugbank_1, 'r') as fp: 
        data = fp.read()

    data2 = data.replace('{http://www.drugbank.ca}','')

    with open(chk_drugbank_2, 'w') as fp: 
        fp.write(data2)
    
    with open(chk_drugbank_2, 'r') as fp: 
        data = json.load(fp)

    df = pd.DataFrame(data['drugbank']['drug'])
    df.to_csv(chk_drugbank_3, sep='\t', index=False)
    df.head(10).to_csv(chk_drugbank_4, sep='\t', index=False)


def load_table():
    process_step1()
    process_step2()

    return pd.read_table(chk_drugbank_3)


def load_json():
    process_step1()
    process_step2()

    with open(chk_drugbank_2, 'r') as fp: 
        data = json.load(fp)


def load_structures():    
    local_filename = download_v2(
        remote_path=url_drugbank_structure,
        download_dir=data_dir
        )

    structure_db = []

    for mol in pybel.readfile('sdf', local_filename):

        dictv = {x: y for x, y 
            in zip(mol.data.keys(), mol.data.values())
            }

        structure_db.append( dictv ) 

    return structure_db





