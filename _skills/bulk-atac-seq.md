---
title: "Bulk ATAC-seq"
collection: skills
permalink: /skills/bulk-atac-seq
excerpt: 'Bulk ATAC-seq for ILCs. <br/><img src="/images/ILC-atac.png">'
date: 2023-11-24
---

Background
======
## Bulk-ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing)
Bulk-ATAC-seq is a powerful high-throughput sequencing technique used to profile chromatin accessibility across the genome. This method leverages a hyperactive Tn5 transposase to insert sequencing adapters into open regions of chromatin, enabling the identification of regulatory elements such as promoters, enhancers, and other DNA-protein interaction sites. As a skilled practitioner of Bulk-ATAC-seq, I am proficient in experimental design, library preparation, and downstream bioinformatics analysis, including peak calling, motif enrichment analysis, and integration with other omics datasets (e.g., RNA-seq, ChIP-seq). My expertise in this technique has enabled me to investigate epigenetic mechanisms underlying gene regulation in various biological contexts, contributing to a deeper understanding of cellular differentiation, disease states, and transcriptional networks. By combining Bulk-ATAC-seq with advanced computational tools, I aim to uncover novel insights into the dynamic interplay between chromatin architecture and gene expression.<br/>

A significant application of my work involves leveraging Bulk-ATAC-seq to study the epigenetic mechanisms governing immune cell differentiation and function. Specifically, I have employed this technique to investigate the chromatin dynamics of T helper (Th) cells and Innate Lymphoid Cells (ILCs), two critical components of the immune system. In Th cells, Bulk-ATAC-seq has allowed me to dissect the regulatory networks underlying their subset differentiation (e.g., Th1, Th2, Th17, and Treg cells) and their roles in immune responses and diseases. Similarly, in ILCs, I have utilized this approach to explore the epigenetic regulation of their development and functional plasticity, shedding light on their contributions to innate immunity and tissue homeostasis.<br/>

Result
======
The assay for transposase-accessible chromatin with sequencing (ATAC-Seq) is a popular method for determining chromatin accessibility across the genome. By sequencing regions of open chromatin, ATAC-Seq can help you uncover how chromatin packaging and other factors affect gene expression. ATAC-Seq does not require prior knowledge of regulatory elements, making it a powerful epigenetic discovery tool. It has been used to better understand chromatin accessibility, transcription factor binding, and gene regulation in complex diseases, embryonic development, T-cell activation, and cancer. ATAC-Seq can be performed on bulk cell populations or on single cells at high resolution.<br/>

**ILC Il18r ATAC-seq IGV GSE116093**<br/><img src="/images/ILC-atac.png"><br/>
**QC**<br/><img src="/images/atac3.png"><br/>
**Other ATAC**<br/><img src="/images/atac4.png"><br/>



Method
======

**Schematic overview of ATAC-seq protocol**<img src="/images/atac1.png"><br/>
**Overview of the steps of ATAC-seq data analysis**<img src="/images/atac2.png"><br/>


## Preparation
Download the package of atac-seq analysis for Linux.(FastQC, Bowtie2, Samtools, MACS2, IGV, picard, et al.)

### Step 1. QC and Mapping
Standard ATAC-seq pipelines will take the FASTQ files as input and perform a series of QC and data cleaning steps followed by alignment to a reference genome. Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data. Bowtie 2 is geared toward aligning relatively short sequencing reads to long genomes. (Based trimmomatic/bowtie2)
```Linux
cd './fastq'
ls *_R1.fq.gz | while read id
do 
	headname=${id%_R1.fq.gz} 
	trimmomatic PE -threads 40 ${headname}_R1.fq.gz ${headname}_R2.fq.gz ${headname}.clean_1.fastq ${headname}.clean.unpaired_1.fastq ${headname}.clean_2.fastq ${headname}.clean.unpaired_2.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
done

if [ ! -d ../mapping ]
then
	mkdir ../mapping
fi


ls *.clean_1.fastq |while read id
do
    headname=${id%.clean_1.fastq} 
    bowtie2 --threads 40 \
	-x /data/labShare/Sequencing/GENOME/mm10/mm10 \
	-1 ${headname}.clean_1.fastq \
	-2 ${headname}.clean_2.fastq \
	-S ../mapping/${headname}.sam
	rm ${headname}.clean_1.fastq
	rm ${headname}.clean_2.fastq
	rm ${headname}.clean.unpaired_1.fastq
	rm ${headname}.clean.unpaired_2.fastq
done
```
**fastQC**<img src="/images/atac3.png"><br/>
 
### Step 2. SamToBam & filter low-quality mapping
Samtools is a suite of programs for interacting with high-throughput sequencing data.
```Linux
cd '../mapping'
ls *.sam |while read id
do
    samtools view -@ 40 -b -F 4 -S ${id} -o ./${id%.sam}.bam
    samtools sort -@ 40 ./${id%.sam}.bam -o ./${id%.sam}.sorted.bam
	rm ${id}
done
```

### Step 3. de duplicates reads
Picard (MarkDuplicates) and SAMTools (rmdup) are the two main softwares used for PCR duplicate removal.
```Linux
if [ ! -d ../dedup ]
then
	mkdir ../dedup
fi


ls *.sorted.bam |while read id
do
	picard MarkDuplicates REMOVE_DUPLICATES=true I=./${id} O=../dedup/${id%.sorted.bam}.dedup.bam M=../dedup/${id%.sorted.bam}.marked_dup_metrics.txt
done
```

### Step 4. Convert to BigWig and Calling Peaks
```Linux
cd '../dedup'
if [ ! -d ../bg ]
then
	mkdir ../bg
fi
ls *.dedup.bam |while read id
do
    samtools index ./${id}
    bamCoverage -b ./${id} -o ../bg/${id%.dedup.bam}.bw -bs 1 --normalizeUsing BPM -p 40
done

if [ ! -d ../macs2 ]
then
	mkdir ../macs2
fi

ls *.dedup.bam |while read id
do
	macs2 callpeak -t ${id} -n sample --shift -75 --extsize 150 --nomodel -B --SPMR -g mm --outdir ../macs2/${id%.dedup.bam}.Macs2_out -q 0.01
	
done
```

### Step 5. Motif analysis
Peak calling is a computational method used to identify areas in a genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing or MeDIP-seq experiment. A commonly used tool for identifying transcription factor binding sites is named Model-based Analysis of ChIP-seq (MACS). The MACS algorithm captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. Although it was developed for the detection of transcription factor binding sites it is also suited for larger regions.
```Linux
if [ ! -d ../peak_motif ]
then
mkdir ../peak_motif
fi

ls ../macs2/*Macs2_out/sample_peaks.narrowPeak | while read id
do
findMotifsGenome.pl ${id} mm10 ../peak_motif/ -size given -len 8,10,12 -p 40
done
```

### One Step code
(1)Modify the Fastq file storage location; (2)Pay attention to the species of the reference. (3) Using sh for the one-step analysis
```Linux
#!/bin/bash


# 1.mapping
cd './fastq'
ls *_R1.fq.gz | while read id
do 
	headname=${id%_R1.fq.gz} 
	trimmomatic PE -threads 40 ${headname}_R1.fq.gz ${headname}_R2.fq.gz ${headname}.clean_1.fastq ${headname}.clean.unpaired_1.fastq ${headname}.clean_2.fastq ${headname}.clean.unpaired_2.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
done

if [ ! -d ../mapping ]
then
	mkdir ../mapping
fi


ls *.clean_1.fastq |while read id
do
    headname=${id%.clean_1.fastq} 
    bowtie2 --threads 40 \
	-x /data/labShare/Sequencing/GENOME/mm10/mm10 \
	-1 ${headname}.clean_1.fastq \
	-2 ${headname}.clean_2.fastq \
	-S ../mapping/${headname}.sam
	rm ${headname}.clean_1.fastq
	rm ${headname}.clean_2.fastq
	rm ${headname}.clean.unpaired_1.fastq
	rm ${headname}.clean.unpaired_2.fastq
done



# 2.SamToBam & filter low-quality mapping
cd '../mapping'
ls *.sam |while read id
do
    samtools view -@ 40 -b -F 4 -S ${id} -o ./${id%.sam}.bam
    samtools sort -@ 40 ./${id%.sam}.bam -o ./${id%.sam}.sorted.bam
	rm ${id}
done

# 3.de duplicates reads
if [ ! -d ../dedup ]
then
	mkdir ../dedup
fi


ls *.sorted.bam |while read id
do
	picard MarkDuplicates REMOVE_DUPLICATES=true I=./${id} O=../dedup/${id%.sorted.bam}.dedup.bam M=../dedup/${id%.sorted.bam}.marked_dup_metrics.txt
done

# 4.convert to BigWig
cd '../dedup'
if [ ! -d ../bg ]
then
	mkdir ../bg
fi
ls *.dedup.bam |while read id
do
    samtools index ./${id}
    bamCoverage -b ./${id} -o ../bg/${id%.dedup.bam}.bw -bs 1 --normalizeUsing BPM -p 40
done

# 5.calling peaks
if [ ! -d ../macs2 ]
then
	mkdir ../macs2
fi

ls *.dedup.bam |while read id
do
	macs2 callpeak -t ${id} -n sample --shift -75 --extsize 150 --nomodel -B --SPMR -g mm --outdir ../macs2/${id%.dedup.bam}.Macs2_out -q 0.01
	
done

# 6.Motif analysis
if [ ! -d ../peak_motif ]
then
mkdir ../peak_motif
fi

ls ../macs2/*Macs2_out/sample_peaks.narrowPeak | while read id
do
findMotifsGenome.pl ${id} mm10 ../peak_motif/ -size given -len 8,10,12 -p 40
done


echo ""
echo "WORK DONE!"
```
