---
title: "Bulk RNA-seq"
collection: skills
permalink: /skills/bulk-rna-seq
excerpt: 'Bulk RNA-seq for DSS_AOMDSS (GSE57569, just for code)'
date: 2024-01-24
---

Result
======

**Bulk RNA sequencing workflow**<img src="/images/Skills-Bulk-rna-seq/upstream-workflow.png"><br/>
Bulk RNA sequencing is the method of choice for transcriptomic analysis of pooled cell populations, tissue sections, or biopsies. It measures the average expression level of individual genes across hundreds to millions of input cells and is useful to get a global idea of gene expression differences between samples.<br/>

**Upstream_Result1**<img src="/images/Skills-Bulk-rna-seq/upstream1.png"><br/>

**Upstream_Result2**<img src="/images/Skills-Bulk-rna-seq/upstream2.png"><br/>
Counting result is most important, using this by R to get the result (Heatmap, Volcano Plot, GO/KEGG, GSEA, and other interested downstream figures)<br/>

**Downstream_Result (for example)**<img src="/images/Skills-Bulk-rna-seq/downstream.png"><br/>

**Downstream_Result (GSE19188 PCA)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-1.png"><br/>

**Downstream_Result (GSE19188 UMAP)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-2.png"><br/>

**Downstream_Result (GSE19188 Volcano Plot)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-3.png"><br/>

**Downstream_Result (GSE19188 Heat Map)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-4.png"><br/>

**Downstream_Result (GSE19188 GO)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-5.png"><br/>

**Downstream_Result (GSE19188 KEGG)**<img src="/images/Skills-Bulk-rna-seq/GSE19188-6.png"><br/>



Method
======
### Code
```Linux
#!/bin/bash

# 2.mapping
cd './fastq'
if [ ! -d ../mapping ]
then
	mkdir ../mapping
fi

ls *_1.fastq.gz |while read id
do
    headname=${id%_1.fastq.gz} 
    hisat2 -p 70 --dta -x /data/labShare/Sequencing/GENOME/hisat2/grcm38_tran/genome_tran -1 ${headname}_1.fastq.gz -2 ${headname}_2.fastq.gz -S ../mapping/${headname}.sam
done



# 3.SamToBam & filter low-quality mapping
cd '../mapping'
ls *sam |while read id
do
    samtools view -@ 70 -b -F 4 -S ${id} -o ./${id%.sam}.bam
    samtools sort -@ 70 ./${id%.sam}.bam -o ./${id%.sam}.sorted.bam
	rm ${id}
	samtools index ./${id%.sam}.sorted.bam
done

# 4.counting
if [ ! -d ../count ]
then
	mkdir ../count
fi
```

 
