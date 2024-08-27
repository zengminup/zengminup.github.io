---
title: "Bulk ATAC-seq"
collection: skills
permalink: /skills/bulk-atac-seq
excerpt: 'Bulk ATAC-seq for ILCs'
date: 2023-11-24
---

Result
======
**ILC Il18r ATAC-seq IGV GSE116093**<img src="/images/ILC-atac.png"><br/>

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
