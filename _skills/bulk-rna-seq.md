---
title: "Bulk RNA-seq"
collection: skills
permalink: /skills/bulk-rna-seq
excerpt: 'Bulk RNA-seq for DSS_AOMDSS (GSE57569, just for code)<br/>
<img src="/images/Skills-Bulk-rna-seq/upstream-workflow.png" style="width:250px; height:600px;" align="center">'
date: 2024-01-24
---

Background
======
## Bulk RNA-Seq (Bulk RNA Sequencing)
Bulk RNA-Seq is a powerful high-throughput sequencing technology that enables comprehensive profiling of the transcriptome, providing insights into gene expression patterns across entire populations of cells. By capturing and sequencing RNA molecules from a tissue or cell population, this technique quantifies transcript abundance, identifies differentially expressed genes, and reveals splicing variants, thereby offering a holistic view of cellular states and functions.<br/>

The core principle of Bulk RNA-Seq involves the extraction of total RNA, conversion to cDNA, and subsequent sequencing using next-generation sequencing (NGS) platforms. Bioinformatics pipelines are then employed to align reads, quantify gene expression, and perform downstream analyses such as pathway enrichment and network modeling. This approach is particularly valuable for addressing complex biological questions, including the identification of biomarkers, understanding disease mechanisms, and uncovering regulatory networks.<br/>

In the field of immunology, Bulk RNA-Seq has emerged as a transformative tool for dissecting immune responses, characterizing immune cell populations, and elucidating the molecular mechanisms underlying immune-related diseases. For instance, it has been instrumental in profiling immune cell activation states, identifying cytokine signatures, and uncovering novel immune checkpoints in conditions such as cancer, autoimmune disorders, and infectious diseases. By integrating Bulk RNA-Seq data with other omics datasets, researchers can gain a deeper understanding of the immune system's complexity and develop targeted therapeutic strategies.<br/>

My expertise in Bulk RNA-Seq encompasses experimental design, library preparation, data analysis, and interpretation, with a particular focus on its application to immunological research. I am proficient in utilizing state-of-the-art computational tools and statistical methods to extract meaningful biological insights from large-scale transcriptomic datasets, contributing to advancements in both basic and translational immunology.<br/>

Result
======

**Bulk RNA Sequencing Workflow**<img src="/images/Skills-Bulk-rna-seq/upstream-workflow.png"><br/>
Bulk RNA sequencing is the method of choice for transcriptomic analysis of pooled cell populations, tissue sections, or biopsies. It measures the average expression level of individual genes across hundreds to millions of input cells and is useful to get a global idea of gene expression differences between samples.<br/>

**Upstream Result1**<img src="/images/Skills-Bulk-rna-seq/upstream1.png"><br/>

**Upstream Result2**<img src="/images/Skills-Bulk-rna-seq/upstream2.png"><br/>
Counting result is most important, using this by R to get the result (Heatmap, Volcano Plot, GO/KEGG, GSEA, and other interested downstream figures)<br/>

**Downstream Result (for example)** <img src="/images/Skills-Bulk-rna-seq/downstream.png"><br/>

**Downstream Result (GSE19188 PCA)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-1.png"><br/>

**Downstream Result (GSE19188 UMAP)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-2.png"><br/>

**Downstream Result (GSE19188 Volcano Plot)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-3.png"><br/>

**Downstream Result (GSE19188 Heat Map)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-4.png"><br/>

**Downstream Result (GSE19188 GO)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-5.png"><br/>

**Downstream Result (GSE19188 KEGG)** <img src="/images/Skills-Bulk-rna-seq/GSE19188-6.png"><br/>



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

 
