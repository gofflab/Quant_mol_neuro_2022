#!/usr/bin/env python

import sys
import pandas as pd
import os
import anndata as ad
import scipy.io as sio
from scipy import sparse
import subprocess as sp


dataDir = 'data/'
GEOAcc = 'GSE125708'

# Download the data from GEO
url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={GEOAcc}&format=file"
sp.call(['wget', url, '-O', dataDir + f'{GEOAcc}_RAW.tar'])

# Extract the tar file
sp.call(['tar', '-xvf', dataDir + f'{GEOAcc}_RAW.tar', '-C', dataDir])

# gunzip the .gz files in dataDir
for file in os.listdir(dataDir):
    if file.endswith('.gz'):
        sp.call(['gunzip', dataDir + file])

# Loop through all the '*.txt' files in dataDir and create a list of annData objects from their contents
adatas = []
samples = []
for file in os.listdir(dataDir):
    if file.endswith('.txt'):
        print('Reading in ' + file)
        #get basename of the file after removing .txt
        basename = os.path.splitext(file)[0]
        samples.append(basename)
        adata = ad.read_csv(dataDir + file, delimiter='\t')
        adata= adata.transpose()
        adata.obs['GSM'] = basename.split("_")[0]
        adata.obs['sample'] = basename.split("_")[1]
        if basename.split("_")[1].startswith("O"):
            adata.obs['age'] = "old"
        else:
            adata.obs['age'] = "young"
        adatas.append(adata)

all = ad.concat(adatas,merge="same")
sparse_X = sparse.csr_matrix(all.X)
all.X = sparse_X

# Write the expression matrix to a matrix market file
sio.mmwrite(f"{dataDir}/GSE129788.mtx",all.X)

# Write the obs to a csv file
all.obs.to_csv(f"{dataDir}/GSE129788_cell_info.csv")

# Write the var to a csv file
all.var.to_csv(f"{dataDir}/GSE129788_gene_info.csv")
