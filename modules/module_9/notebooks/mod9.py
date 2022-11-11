#%%
import polars as pl
from pydantic import BaseModel, root_validator, validator

#%%
GEO = "GSE63482"
data = pl.read_csv(
    "https://github.com/gofflab/Quant_mol_neuro_2022/raw/main/modules/module_9/notebooks/GSE63482_Expression_matrix.tsv",
    sep=" ",
)

gene_metadata = pl.DataFrame({"gene_id": data["gene_id"]})

sp = [x.split("_") for x in data.columns[1:]]
sample_metadata = pl.DataFrame(
    dict(
        age=[x[0] for x in sp], cell_type=[x[1] for x in sp], sampleId=data.columns[1:]
    )
)

# %%
class SummarizedExperiment(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    assays: dict[str, pl.DataFrame]
    col_data: pl.DataFrame
    row_data: pl.DataFrame

    @root_validator(pre=True)
    def check_assays(cls, values):
        assays = values["assays"]
        for v in assays.values():
            if v.shape[0] != values["row_data"].shape[0]:
                raise ValueError(
                    f"Assay has {v.shape[0]} rows, but rowData has {values['row_data'].shape[0]} rows"
                )

            # if v.shape[1] != values["col_data"].shape[0]:
            #     raise ValueError(
            #         f"Assay has {v.shape[1]} columns, but colData has {values['col_data'].shape[0]} rows"
            #     )
        return values


se = SummarizedExperiment(
    assays={"fpkm": data}, col_data=sample_metadata, row_data=gene_metadata
)

genes = ["Sox2", "p53"]
sample = "E15_cpn"
fpkm = se.assays["fpkm"]

fpkm.filter(pl.col("gene_id").str.contains("^Sox2$|^p53$")).select(["gene_id", sample])

# %%
fpkm.filter(pl.col("gene_id").str.contains("^Hox")).select(["gene_id", sample])

# %%
expressed = fpkm.filter(fpkm.select(pl.exclude("gene_id")).sum(axis=1) > 0)
import matplotlib.pyplot as plt
import seaborn as sns

# %%
from sklearn.decomposition import PCA

sns.set()


pca = PCA(n_components=2).fit(expressed.select(pl.exclude("gene_id")).transpose())
pca_df = pl.DataFrame(
    {
        "PC1": pca.components_[0],
        "PC2": pca.components_[1],
        "cell_type": sample_metadata["cell_type"],
    }
)
sns.scatterplot(data=pca_df.to_pandas(), x="PC1", y="PC2", hue="cell_type")

#%%
import seaborn.objects as so

fpkm = se.assays["fpkm"].to_pandas()

# %%
p = so.Plot(fpkm, x="E15_cpn", y="E16_cpn").add(so.Dot())
p
# %%
p.add(so.Line(), so.PolyFit())

# %%
p.scale(x="log", y="log")
# %%

genes = ["Sox2", "p53", "Fezf2"]

gene_data = se.assays["fpkm"].filter(
    pl.col("gene_id").str.contains("|".join([f"^{x}$" for x in genes]))
)

data_melted = gene_data.melt(
    id_vars=["gene_id"], variable_name="sampleId", value_name="fpkm"
)
data_melted = data_melted.join(sample_metadata, on="sampleId")


# %%
sns.violinplot(data=data_melted.to_pandas(), x="cell_type", y="fpkm", hue="cell_type")
# %%
(
    so.Plot(data_melted.to_pandas(), x="age", y="fpkm", color="cell_type")
    .facet(row="gene_id", col="cell_type")
    .add(so.Line())
    .scale(y="log")
)
# %%
