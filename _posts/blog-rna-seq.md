---
title: "Single Cell RNA-seq Tutorial (PBMC3K Example)"
permalink: /posts/blog-rna-seq/
date: 2023-09-01
excerpt: 'Single Cell RNA-seq analysis for mice intestinal lamina propria and Peyer’s patches for CD45+ immune cells (GSE124880, Immunity)'
tags:
  - Bioinformatics
---

Method
======

### Step 0. Import Seurat package
First of all, loading package
```R
library(dplyr)
library(Seurat)
library(patchwork)
```
 
 ### Step 1. Create a Seurat object
 ```R
counts <- Read10X(data.dir = "data/DS1/")
seurat <- CreateSeuratObject(counts, project="DS1")
```