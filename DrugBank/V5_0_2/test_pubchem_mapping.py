import pandas as pd 
import os
from os.path import join,dirname,exists
from ipdb import set_trace
import pubchempy
from tqdm import tqdm 

data_dir = dirname(__file__) + '/data'
os.system('mkdir -p {}'.format(data_dir))

#
# input data 
#

dbid_name_inchi = data_dir+'/chk-dbid-name-inchi.tsv'

#
# output data
#

output_dbid_cid_mapping = data_dir+'/chk-pubchem-mapping.tsv'


def test_main():
    drugbank_df = pd.read_table(dbid_name_inchi)
    rows = list()
    k = 0 
    for idx in tqdm(drugbank_df.index):
        row = drugbank_df.loc[idx]
        row['PubChemCID'] = None

        if pd.isnull(row.InChI):
            rows.append(row) 
            continue

        try:
            # set_trace()
            compounds = pubchempy.get_compounds(row.InChI, namespace='inchi')
            # compounds = pubchempy.get_compounds(row.InChIKey, namespace='inchikey')
        except pubchempy.BadRequestError:
            print('BadRequestError', row)
            rows.append(row)
            continue
        
        try:
            compound, = compounds
        except ValueError:
            print(row, compounds)
            rows.append(row)
            continue

        row['PubChemCID'] = '%d' % (compound.cid)

        rows.append(row)

        # k = k + 1 
        # if k > 3: 
        #     break 
            
    df2 = pd.DataFrame(rows)
    df2.to_csv(output_dbid_cid_mapping, sep='\t', index=False)
