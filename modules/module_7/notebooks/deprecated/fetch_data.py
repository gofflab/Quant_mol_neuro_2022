import gzip
import tarfile
from pathlib import Path

import anndata as ad
import pandas as pd
import requests
import scipy.io as sio
from scipy import sparse

dataDir = Path("data")
dataDir.mkdir(exist_ok=True)
GEOAcc = "GSE125708"

# Download the data from GEO
url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={GEOAcc}&format=file"
if not (dataDir / f"{GEOAcc}.tar").exists():
    print("Downloading data...")
    r = requests.get(url, allow_redirects=True)
    (dataDir / f"{GEOAcc}.tar").write_bytes(r.content)

# Extract the tar file
with tarfile.open(dataDir / f"{GEOAcc}.tar") as tar:
    tar.extractall(dataDir)

# gunzip the .gz files in dataDir
for f in dataDir.glob("*.gz"):
    with gzip.open(f, "rb") as gz:
        (dataDir / f.stem).write_bytes(gz.read())
    f.unlink()

# Loop through all the '*.txt' files in dataDir and create a list of annData objects from their contentssamples = {}
samples = {}
for file in dataDir.glob("*.txt"):
    print("Reading in " + file.stem)
    names = file.stem.split("_")

    # Format data such that gene expression is in columns
    df = pd.read_csv(file, delimiter=" ", header=[0]).T
    df.rename(columns=df.iloc[0], inplace=True)
    df.drop(df.index[0], inplace=True)
    df.index = pd.Index([names[0]], dtype="object")
    df = df.astype("float32")

    adata = ad.AnnData(
        X=sparse.csr_matrix(df.values),
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )

    adata.obs["sample"] = names[1]
    if names[1].startswith("A"):
        adata.obs["tissue"] = "retina"
    elif names[1].startswith("J"):
        adata.obs["tissue"] = "whole_brain"
    else:
        raise ValueError("Unknown tissue type")

    # get basename of the file after removing .txt
    samples[file.stem] = adata

df = ad.concat(samples.values(), merge="same")

# Write the expression matrix to a matrix market file
sio.mmwrite(dataDir / "GSE129788.mtx", df.X)

# Write the obs to a csv file
df.obs.to_csv(dataDir / "GSE129788_cell_info.csv")

# Write the var to a csv file
df.var.to_csv(dataDir / "GSE129788_gene_info.csv")
