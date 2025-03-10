---
title: "Spatial Transcriptomics Analysis Tutorial (Scanpy)"
permalink: /posts/blog-spatial-transcriptomics/
date: 2024-10-01
excerpt: 'Spatial Transcriptomics analysis tutorial. Example using Python and Scanpy. <br/>
<img src="/images/spatial_scanpy/spatial2.png" style="width:400px; height:650px;">'
tags:
  - Bioinformatics
---


Background
======  
Spatial transcriptomics is a powerful technique that enables the study of gene expression patterns in tissues with spatial resolution. It provides insights into how cells interact within their native environments, offering a comprehensive understanding of tissue architecture, cellular localization, and the molecular mechanisms that drive tissue function and disease progression.<br>

The fundamental principle behind spatial transcriptomics is the ability to capture gene expression data while preserving the spatial organization of tissues. This technique uses high-throughput sequencing to measure gene expression directly from tissue sections, allowing researchers to map the transcriptome to its spatial context. Unlike bulk RNA sequencing, which averages gene expression across millions of cells, spatial transcriptomics preserves the tissue's spatial context, enabling the identification of spatially distinct cellular populations and their interactions.<br>

This method is particularly powerful in studying complex tissues, such as the brain, tumors, or developmental stages, where tissue organization plays a critical role in the biological processes. By integrating transcriptomics data with spatial information, it becomes possible to uncover gene expression patterns that are linked to tissue-specific functions, disease progression, and therapeutic responses.<br>

In recent years, advancements in computational tools and experimental protocols have made spatial transcriptomics more accessible. One such tool is Scanpy, a Python-based package that is widely used for the analysis and visualization of single-cell RNA-seq data. Scanpy has expanded its capabilities to accommodate spatial transcriptomics data, offering tools for preprocessing, clustering, visualization, and spatial pattern analysis.<br>

Overall, spatial transcriptomics is transforming the way we study tissues and their molecular profiles. By combining spatial and transcriptomic data, this technique holds the potential to reveal new biological insights, identify spatial biomarkers, and inform the development of novel therapeutic strategies for various diseases.<br>

Results
======
**QC** <img src="/images/spatial_scanpy/QC1.png"><br/>
**UMAP** <img src="/images/spatial_scanpy/UMAP.png"><br/>
**Spatial1** <img src="/images/spatial_scanpy/spatial1.png"><br/>
**Spatial2** <img src="/images/spatial_scanpy/spatial2.png"><br/>
**Spatial3** <img src="/images/spatial_scanpy/spatial3.png"><br/>
**Spatial4** <img src="/images/spatial_scanpy/spatial4.png"><br/>
**Spatial5** <img src="/images/spatial_scanpy/spatial5.png"><br/>
**Spatial6** <img src="/images/spatial_scanpy/spatial6.png"><br/>



Method
======

The analysis of spatial transcriptomics data involves several key steps, from preprocessing raw data to visualizing gene expression patterns in the context of tissue spatial organization. Below is a general method to approach the analysis using Scanpy:

#### 1. **Data Preprocessing**
   - **Input Data:** The raw spatial transcriptomics data typically includes gene expression matrices, spatial coordinates of the tissue sections, and possibly metadata related to the tissue or experimental conditions.
   - **Quality Control:** Perform quality control to remove low-quality spots (e.g., spots with too few genes detected or too many mitochondrial genes).
   - **Normalization:** Normalize the gene expression data to account for differences in sequencing depth across different tissue spots or sections.
   - **Filtering:** Filter out genes that are not expressed in the tissue, as well as genes that are expressed across all tissue spots, to reduce noise and focus on genes that show variable expression.

#### 2. **Dimensionality Reduction**
   - **PCA (Principal Component Analysis):** Reduce the dimensionality of the gene expression matrix to make the data more manageable. PCA is often used to capture the major sources of variation in the data.
   - **Spatially Informed Embedding:** Use spatially aware dimensionality reduction methods, such as UMAP or t-SNE, that take spatial coordinates into account when reducing dimensionality. These methods can help to visualize clusters or patterns in gene expression that are related to spatial location.

#### 3. **Clustering and Cell-Type Identification**
   - **Clustering:** Use unsupervised clustering algorithms, such as Louvain or Leiden clustering, to identify distinct groups of spots with similar gene expression profiles. These clusters may correspond to different tissue regions or cell types.
   - **Cell-Type Annotation:** If available, map cell-type markers to the clusters to identify the cell types represented in each cluster. This can be done by using known marker genes for specific cell types or by performing a supervised learning approach with labeled data.
   
#### 4. **Spatial Pattern Analysis**
   - **Spatially Variable Genes (SVGs):** Identify genes that exhibit spatially variable expression across the tissue, indicating that their expression is influenced by the local tissue environment.
   - **Gene Expression Mapping:** Map the expression of specific genes or gene sets onto the tissue’s spatial coordinates. This provides insight into how gene expression is distributed across the tissue and how it relates to tissue architecture.
   - **Differential Expression:** Perform differential expression analysis to identify genes that are differentially expressed between spatially distinct regions of the tissue.

#### 5. **Visualization**
   - **Spatial Plots:** Generate spatial plots to visualize the expression of genes or clusters across the tissue. Common visualizations include heatmaps, dot plots, or 2D projections that incorporate both gene expression and spatial information.
   - **Interactive Visualization:** Use interactive plotting tools to explore the spatial distribution of gene expression and clusters. This allows for better exploration of the data and identification of patterns that are not immediately obvious.
   
#### 6. **Integration with Other Data**
   - **Integrate with Single-Cell Data:** Combine spatial transcriptomics data with single-cell RNA-seq data to enhance cell-type identification and understand how cell types are distributed across tissue sections.
   - **Cross-Tissue Comparison:** Compare spatial transcriptomics data across different tissue types or conditions to uncover spatially specific gene expression patterns that may be associated with disease or developmental stages.

#### 7. **Interpretation and Biological Insights**
   - **Gene Function Analysis:** Investigate the biological function of genes that are spatially variable or differentially expressed across the tissue. This can provide insights into the molecular mechanisms underlying tissue organization, disease processes, or therapeutic responses.
   - **Pathway Analysis:** Perform pathway enrichment analysis on the identified genes to gain a deeper understanding of the biological processes that are spatially regulated in the tissue.

By following these steps, spatial transcriptomics data can be analyzed and visualized in a way that reveals valuable insights into the spatial organization of gene expression and its relationship to tissue function and pathology.


```Linux
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Print the version of Scanpy to ensure the correct environment configuration
sc.logging.print_versions()

# Set the figure parameters, including the background color and figure size
sc.set_figure_params(facecolor="white", figsize=(8, 8))

# Set the verbosity level for Scanpy logs to 3 for detailed output
sc.settings.verbosity = 3

# Load the Visium SGE human lymph node dataset as an example, containing spatial transcriptomics data
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")

# Ensure gene names are unique to avoid conflicts with duplicated gene names
adata.var_names_make_unique()

# Annotate mitochondrial genes by checking whether the gene names start with "MT-"
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate quality control metrics, including the percentage of mitochondrial gene counts for each cell
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Create subplots to visualize the distribution of QC metrics
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
# Plot the distribution of total counts (total gene expression) across all cells
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
# Plot the distribution of total counts, excluding extreme values (<10,000)
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
# Plot the distribution of the number of genes detected in each cell
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
# Plot the distribution of genes detected in each cell, excluding extreme values (<4000)
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])

# Filter out cells with low total counts (<5000) and cells with extremely high total counts (>35000)
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)

# Further filter cells with a mitochondrial gene proportion greater than 20%
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
print(f"#cells after MT filter: {adata.n_obs}")

# Filter out genes that are expressed in fewer than 10 cells
sc.pp.filter_genes(adata, min_cells=10)

# Normalize the total counts per cell to account for differences in sequencing depth
sc.pp.normalize_total(adata, inplace=True)

# Log-transform the gene expression data to stabilize variance and facilitate downstream analysis
sc.pp.log1p(adata)

# Identify highly variable genes using the Seurat method, selecting the top 2000 genes
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Perform Principal Component Analysis (PCA) to reduce the dimensionality of the data
sc.pp.pca(adata)

# Compute the neighborhood graph to enable clustering and dimensionality reduction
sc.pp.neighbors(adata)

# Perform Uniform Manifold Approximation and Projection (UMAP) for non-linear dimensionality reduction
sc.tl.umap(adata)

# Perform Leiden clustering to identify cell groups based on gene expression profiles
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)

# Adjust figure size and plot the UMAP embedding, showing total counts, number of genes, and clusters
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

# Adjust figure size and plot the spatial distribution of total counts and gene counts across the tissue section
plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])

# Plot the spatial distribution of clusters across the tissue section
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)

# Crop the image to focus on a specific region of interest and visualize the clusters
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)

# Perform differential gene expression analysis using a t-test to identify marker genes for each cluster
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")

# Plot a heatmap of the top 10 marker genes for cluster 9
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")

# Visualize the spatial distribution of gene expression for the CR2 gene
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])

# Visualize the spatial distribution of genes COL1A2 and SYPL1 with transparency (alpha=0.7)
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)

# Visualize the spatial distribution of genes IL18 and IL18R1 with transparency (alpha=0.7)
sc.pl.spatial(adata, img_key="hires", color=["IL18", "IL18R1"], alpha=0.7)
```
