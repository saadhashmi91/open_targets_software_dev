#!/usr/bin/env python
# coding: utf-8

# In[17]:


import sys
import subprocess

# implement pip as a subprocess:
packages = ['pandas','numpy','pyarrow']
for p in packages:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install',p])


# In[18]:


import os
import pandas as pd
import numpy as np
import io
import fnmatch
import tarfile
from ftplib import FTP,error_perm
from IPython.display import display


# In[3]:


paths_to_df = [
 ('targets.feather','targets.tar.gz','pub/databases/opentargets/platform/21.11/output/etl/parquet/targets')
,('diseases.feather','diseases.tar.gz','pub/databases/opentargets/platform/21.11/output/etl/parquet/diseases')
,('evidence.feather','evidence.tar.gz','pub/databases/opentargets/platform/21.11/output/etl/parquet/evidence/sourceId=eva')
]


# In[4]:


def load_datasets_with_fallback():
    try:
        for feather_path,parquet_path,ftp_path in paths_to_df:
            if os.path.isfile(feather_path):
                yield pd.read_feather(feather_path)
            elif os.path.isfile(parquet_path):
                with tarfile.open(name, "r:*") as tar:
                    df = pd.concat([pd.read_parquet(tar.extractfile(file)) for file in tar.getnames()])
                    df.reset_index().to_feather(feather_path,compression = 'lz4')
                    yield pd.read_feather(feather_path)
            else:
                ftp = FTP('ftp.ebi.ac.uk')
                ftp.login()
                ftp.cwd(path)
                files = [filename for filename in ftp.nlst() if fnmatch.fnmatch(filename, '*.parquet')]
                with tarfile.open(name, "w:gz") as fout:
                    for file_ in files:
                        buffer = io.BytesIO()
                        ftp.retrbinary('RETR ' + str(file_), buffer.write)
                        size = buffer.getbuffer().nbytes
                        print(f"Size of {file_} is: {size}")
                        tf = tarfile.TarInfo(file_)
                        tf.size = size
                        buffer.seek(0)
                        fout.addfile(tf,buffer)
                with tarfile.open(name, "r:*") as tar:
                    df = pd.concat([pd.read_parquet(tar.extractfile(file)) for file in tar.getnames()])
                    df.reset_index().to_feather('diseases.feather',compression = 'lz4')
                    yield pd.read_feather(feather_path)
    except error_perm as resp:
        if str(resp) == "550 No files found":
            print("No files in this directory")
        else:
            raise


# In[5]:


targets_df,diseases_df,evidence_df = list(load_datasets_with_fallback())
display(targets_df)
display(diseases_df)
display(evidence_df)


# In[6]:


evidence_df = evidence_df[['targetId','diseaseId','score']]
c = ['targetId', 'diseaseId']
# trgt_disease_with_median = evidence_df.groupby(c)['score'].median()
trgt_disease_with_medain_top3   = evidence_df.groupby(c,as_index=False).agg(median = pd.NamedAgg(column='score', aggfunc='median'),
                                                      top3   = pd.NamedAgg(column='score', aggfunc= lambda x: (x.nlargest(3).tolist())
                                                                              ))


# In[7]:


joined_targets_df = trgt_disease_with_medain_top3.merge(targets_df[['id','approvedSymbol']],left_on='targetId',right_on='id')

joined_targets_df.drop(['id'],axis=1,inplace=True)
joined_targets_df


# In[8]:


joined_target_diesease_df = joined_targets_df.merge(diseases_df[['id','name']],left_on = 'diseaseId',right_on = 'id')
joined_target_diesease_df.drop(['id'],axis=1,inplace=True)
joined_target_diesease_df


# In[9]:


joined_target_diesease_df.sort_values(by=['median']).to_json('output.json',orient= 'records')

print('Result for part one is now stored as output.json')


# In[10]:


target_disease_df = joined_target_diesease_df[['targetId','diseaseId']]


# In[11]:


result = pd.merge(target_disease_df,target_disease_df, on='diseaseId',suffixes =['','_right'])

result = result[result.targetId != result.targetId_right]

c = ['targetId','targetId_right']

target_target_disease_count = result.groupby(c,as_index=False).count()
target_target_disease_count = target_target_disease_count[target_target_disease_count.diseaseId >= 2]

print('Number of target-target pairs which share a connection to at least two diseases:',len(target_target_disease_count.index))

