---
title: "Single Cell RNA-seq "
collection: skills
permalink: /skills/sc-rna-seq-scanpy
excerpt: 'Single Cell RNA-seq analysis for PBMC3k based Python-Scanpy'
date: 2022-12-15
---

Process
======

### Step 0. Get the sequence data
```Linux
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

 
### Step 1. Import module and sequence data (Using Linux insstead of Windows based Python because of based-environment)
 ```Python
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
results_file = "write/pbmc3k.h5ad"  # the file that will store the analysis results
```

### Step 2. Read matrix into an AnnData object
```Python
 adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
```

### Step 3. Preprocessing and quality control
Show those genes that yield the highest fraction of counts in each single cell, across all cells.
```Python
sc.pl.highest_expr_genes(adata, n_top=20)
```
**Top 20 gene** <img src="/images/scanpy1.png"><br/>

A violin plot of some of the computed quality measures: the number of genes expressed in the count matrix; the total counts per cell; the percentage of counts in mitochondrial genes
```Python
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
```
**Violin QC Plot** <img src="/images/scanpy2.png"><br/>
**QC Plot** <img src="/images/scanpy3.png"><br/>

```Python
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
```
**Highly variable genes** <img src="/images/scanpy5.png"><br/>

Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
```Python
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
```

### Step 4. Principal component analysis
Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
```Python
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca(adata, color="CST3")
```
**PCA Plot** <img src="/images/scanpy6.png"><br/>
Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function sc.tl.louvain() or tSNE sc.tl.tsne(). In our experience, often a rough estimate of the number of PCs does fine.
```Python
sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)
```
**Variance Ratio Plot** <img src="/images/scanpy7.png"><br/>

### Step 5. Computing the neighborhood graph
Let us compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
```Python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["CST3", "NKG7", "PPBP"])
```
**UMAP Plot** <img src="/images/scanpy8.png"><br/>

### Step 6. Clustering the neighborhood graph
```Python
sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
sc.pl.umap(adata, color=["leiden", "CST3", "NKG7"])
adata.write(results_file)
```
**UMAP Plot2** <img src="/images/scanpy9.png"><br/>

### Step 6. Finding marker genes and annotation
Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.
```Python
sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.settings.verbosity = 2  # reduce the verbosity
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(results_file)
sc.tl.rank_genes_groups(adata, "leiden", method="logreg", max_iter=1000)
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
marker_genes = [
    *["IL7R", "CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "CD14"],
    *["LGALS3", "S100A8", "GNLY", "NKG7", "KLRB1"],
    *["FCGR3A", "MS4A7", "FCER1A", "CST3", "PPBP"],
]
adata = sc.read(results_file)
pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
pd.DataFrame(
    {
        group + "_" + key[:1]: result[key][group]
        for group in groups
        for key in ["names", "pvals"]
    }
).head(5)
sc.tl.rank_genes_groups(adata, "leiden", groups=["0"], reference="1", method="wilcoxon")
sc.pl.rank_genes_groups(adata, groups=["0"], n_genes=20)
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)
adata = sc.read(results_file)
sc.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)
sc.pl.violin(adata, ["CST3", "NKG7", "PPBP"], groupby="leiden")
new_cluster_names = [
    "CD4 T",
    "B",
    "FCGR3A+ Monocytes",
    "NK",
    "CD8 T",
    "CD14+ Monocytes",
    "Dendritic",
    "Megakaryocytes",
]
adata.rename_categories("leiden", new_cluster_names)
sc.pl.umap(
    adata, color="leiden", legend_loc="on data", title="", frameon=False, save=".pdf"
)
sc.pl.dotplot(adata, marker_genes, groupby="leiden");
```
**Annotation Plot** <img src="/images/scanpy10.png"><br/>
**Dot Plot2** <img src="/images/scanpy11.png"><br/>
