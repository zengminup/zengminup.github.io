---
title: "Single Cell RNA-seq Tutorial (PBMC3K Example)"
date: 2023-09-01
permalink: /posts/blog-rna-seq/
excerpt: 'Single Cell RNA-seq analysis tutorial for PMBC3K example'
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