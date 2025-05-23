---
title: "Single Cell RNA-seq (Based on R-Seurat)"
collection: skills
permalink: /skills/sc-rna-seq
excerpt: 'Single Cell RNA-seq analysis for mice intestinal lamina propria and Peyer’s patches for CD45+ immune cells (GSE124880, Immunity). <br/>
<img src="/images/GSE124880UMAP.png" style="width:800px; height:550px;" align="center">'
date: 2023-10-12
---

Background
======
## Single Cell RNA-seq (scRNA-seq)
Single Cell RNA-seq is a transformative technology that enables the profiling of gene expression at single-cell resolution, offering unprecedented insights into cellular heterogeneity, developmental trajectories, and functional states within complex tissues. As a proficient user of scRNA-seq, I specialize in experimental design, single-cell library preparation, and advanced computational analysis to dissect the transcriptional landscapes of diverse cell populations. My work extensively utilizes R-Seurat, a powerful and widely adopted toolkit for scRNA-seq data analysis, to perform quality control, dimensionality reduction, clustering, and differential gene expression analysis.<br>

In the field of immunology, I have applied scRNA-seq to investigate the functional diversity and regulatory mechanisms of immune cells, with a particular focus on T helper (Th) cells and Innate Lymphoid Cells (ILCs). By integrating scRNA-seq with other omics approaches, such as Bulk-ATAC-seq, I have uncovered novel insights into the transcriptional programs driving Th cell subset differentiation (e.g., Th1, Th2, Th17, and Treg cells) and the plasticity of ILCs in response to environmental cues. These findings have deepened our understanding of immune cell dynamics in health and disease, including autoimmune disorders, infections, and cancer.<br>

Through the application of scRNA-seq and tools like Seurat, I aim to unravel the complexity of immune cell interactions and signaling networks, paving the way for the development of targeted immunotherapies and precision medicine approaches. My expertise in this area bridges experimental and computational biology, enabling the discovery of novel biomarkers and therapeutic targets in immunology and beyond.<br>

Result
======
**UMAP Figure**<br><img src='/images/GSE124880UMAP.png'><br>

**Cell Annotation Figure**<br><img src='/images/ILC.png'><br><img src='/images/T.png'>
 
**VlnPlot Figure**<br><img src='/images/GSE124880VlnPlot.png'><br>
 
**FeaturePlot Figure**<br><img src='/images/GSE124880FeaturePlot.png'><br>
 
**DotPlot Figure**<br><img src='/images/GSE124880DotPlot.png'><br>


Method
======
R package and code can be shown in github and Blogs. <br/>
[Github-Link](https://github.com/zengminup/SC-RNA-seq)<br/>
[Blog-Link](https://zengminup.github.io/posts/blog-rna-seq/)<br/>
<br/>
## Preparation
This tutorial assumes that the sequencing data preprocessing steps, including base calling, mapping and read counting, have been done. 10x Genomics has its own analysis pipeline [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) for data generated with the 10x Genomics Chromium Single Cell Gene Expression Solution. At the end of the Cell Ranger pipeline, a count matrix is generated. If your scRNA-seq data is generated using another technology (e.g. well-based experiments using Smart-Seq2 and others), the Cell Ranger pipeline is likely unapplicable, and you will have to find another solution to generate the count matrix.
```Linux
mkdir ~/yard/run_cellranger_count
cd ~/yard/run_cellranger_count
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf pbmc_1k_v3_fastqs.tar
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
cellranger count --id=run_count_1kpbmcs \
   --fastqs=/data/zengmin/github/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/data/zengmin/github/refdata-gex-GRCh38-2020-A
```

### Step 0. Import Seurat package
First of all, loading Seurat package
```R
library(dplyr)
library(Seurat)
library(patchwork)
```
 
 ### Step 1. Create a Seurat object
 ```R
seurat.data <- Read10X(data.dir = "/data/zengmin/github/run_count_1kpbmcs/outs/raw_feature_bc_matrix/")
seurat <- CreateSeuratObject(counts = seurat.data, project = "seurat3k", min.cells = 3, min.features = 200)
```

### Step 2. Quality control
After creating the Seurat object, the next step is to do quality control on the data. The most common quality control is to filter out
1. Cells with too few genes detected.
2. Cells with too many genes detected.
3. Cells with high mitochondrial transcript percentage.
```R
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

<img src="/images/QC.png"><br/>


And as one would expect, number of detected genes and number of detected transcripts are well correlated across cells while mitochondrial transcript percentage is not.
```R
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

<img src="/images/QC2.png"><br/>

For instance, for this data set, a detected gene number between 500 and 5000, and a mitochondrial transcript percentage lower than 5% would be quite reasonable, but it is fine to use different thresholds.
```R
seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 5)
```
It is worth to mention that sometimes more QC may need to be applied. DoubletFinder is need for QC, [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder).

### Step 3. Normalization
A normalization step, aiming to make gene expression levels between different cells comparable, is therefore necessary.
```R
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
```

### Step 4. Feature selection for following heterogeneity analysis
In Seurat, or more general in scRNA-seq data analysis, this step usually refers to the identification of highly variable features/genes, which are genes with the most varied expression levels across cells. One can change the number of highly variable features easily by giving the nfeatures option.
One can visualize the result in a variable feature plot, but this is optional.
```R
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat), 10)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```
<img src="/images/Findvar.png"><br/>

### Step 5. Data scaling
Since different genes have different base expression levels and distributions, the contribution of each gene to the analysis is different if no data transformation is performed. This is something we do not want as we don't want our analysis to only depend on genes that are highly expressed. Therefore a scaling is applied to the data using the selected features, just like one usually does in any data science field.
```R
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
```

### (Optional and advanced) Alternative step 3-5: using SCTransform
One problem of doing the typical log-normalization is that is [introduces the zero-inflation artifact](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) into the scRNA-seq data. To better resolve this issue, Hafemeister and Satija introduced an R package ```sctransform```, which uses a regularized negative binomial regression model to normalize scRNA-seq data. Seurat has a wrapper function ```SCTransform```.
```R
seurat <- SCTransform(seurat, variable.features.n = 3000)
```
The ```variable.features.n``` controls the number of highly variable features to identify. One can also add information about which unwanted sources of variation to regress out. For instance,
```R
seurat <- SCTransform(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"), variable.features.n = 3000)
```
This operation combines normalization, scaling and highly variable feature identification so it essentially replaces steps 3-5 from above. Drawbacks of running ```SCTransform``` include
1. It is relatively slow.
2. It makes the normalized expression measurements data-dependent. In the standard procedure, the normalization only relies on the cell itself; in ```SCTransform```, however, information from the other cells in the same data set is involved during normalization. This potentially introduces problems when multiple data sets need to be compared, since the normalized expression measurements of two data sets each individually normalized using ```SCTransform``` are not comparable.
3. There are steps in ```SCTransform``` which involve random sampling to speed up the computation. That means that there is some stochasticity in ```SCTransform``` and the result is slightly different every time, even if it is applied to the same data set.

### Step 6. Linear dimensionality reduction using principal component analysis (PCA)
Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
```R
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
DimPlot(seurat, reduction = "pca") + NoLegend()
```
<img src="/images/pca.png"><br/>

Apart from making an unbiased decision, one can also check which genes are mostly contributing to each of the top PCs . This can be informative if one knows the genes and the biology of the analyzed sample. It provides the opportunity to understand the biological implication of each of the top PCs, so that one can pick those representing useful information.
```R
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
```
<img src="/images/heat3.png"><br/>

### Step 7. Non-linear dimension reduction for visualization
The most commonly used non-linear dimension reduction methods in scRNA-seq data analysis are t-distributed Stochastic Neighbor Embedding (t-SNE) and Uniform Manifold Approximation and Projection (UMAP). Both methods try to place every sample in a low-dimensional space (2D/3D), so that distances or neighborhood relationships between different samples (here cells) in the original space are largely retained in the low-dimensional space. The detailed mathematically descriptions of the two methods are out of the scope of this tutorial, but for those who are interested in, you may check this [video](https://www.youtube.com/watch?v=RJVL80Gg3lA&list=UUtXKDgv1AVoG88PLl8nGXmw) for tSNE, and this [blog](https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668) of Nikolay Oskolkov for UMAP. There are also more methods to create other low-dimensional embeddings for visualization, including but not limiting to [SPRING](https://kleintools.hms.harvard.edu/tools/spring.html), [PHATE](https://phate.readthedocs.io/en/stable/). Now let's focus on tSNE and UMAP which Seurat has included. The top PCs in the PCA analysis are used as the input to create a tSNE and UMAP embedding of the data.
```R
ElbowPlot(seurat)
seurat <- RunUMAP(seurat, dims = 1:10)
seurat <- RunTSNE(seurat, dims = 1:10)
DimPlot(seurat, reduction = "umap")
```
<img src="/images/umap.png" align="centre" /><br/>

### Step 8. Cluster the cells
```R
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.5)
pbmc.markers <- FindAllMarkers(seurat, only.pos = TRUE)
```

### Step 9. Annotate cell clusters
1. Check the expression of canonical cell type and cell state markers in these clusters;
2. Identify signature genes, or marker genes, of each identified cell cluster. Based on the identified cluster marker genes, one can do literature search, enrichment analysis or do experiment (or ask people around) for annotation;
3. For each cluster, compare its gene expression profile with existing reference data.
```R
FeaturePlot(seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK")
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 1, label.size = 4)
```
<img src="/images/feature.png" align="centre" /><br/>
<img src="/images/umap2.png" align="centre" /><br/>

All Code
```R
library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "/brahms/mollag/practice/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```
