all but RNA-seq_4patients.pbs and RNA-seq_4patients_diff.R scripts can be independently run. But all the scripts should run before running other scripts for figures (in figure folders)
First run RNA-seq_4patients.pbs and then RNA-seq_4patients_diff.R should run Other scripts are independent from each other

##########################################################################################
##########################################################################################
##########################################################################################
description about 450k_celllines.rmd

This script finds differentially methylated CpGs from 450k arrays. raw idat Files are in 450k_methylation folder
files related to patients are in patient folder
files related to 2 cell lines are in cellline folder

Reference genome for 450k files is hg19

run the scripts and then files will be generated

It will output differentially methylated CpGs information for other downstream analysis
all the files are in hg19 

##########################################################################################
##########################################################################################
##########################################################################################

description about  RNA-seq_all_patients.Rmd

RNA-seq differential expression of all VHL mutated patients

Differential gene expression between tumor and normal from 29 patients
inputs:
##########################################################################################
rawCountMatrix.tsv

This file contains the read counts of genes for each patient as a matrix (obtained from Scelo, G. et al. Variation in genomic landscape of clear cell renal cell carcinoma across Europe. Nature Communications 5, 5135 (2014).)

first column is the gene name
and other columns are read counts for each patients  for both tumor and normal samples


##########################################################################################
seelcted.patient.data.sample.sheet.xlsx (obtained from Scelo, G. et al. Variation in genomic landscape of clear cell renal cell carcinoma across Europe. Nature Communications 5, 5135 (2014))

This file contains the clinical information about patients

##########################################################################################
output:

ordered down regulated genes paired conventional patients.csv (GRCh 38)

File contains differential gene expression information of genes

Columns

"baseMean expression values","log2FoldChange","Standard Error of log2 fold change ","stat","pvalue","padj"
Read more about DESeq2 output
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

##########################################################################################
##########################################################################################
##########################################################################################

Description about   RNA-seq_4patients.pbs

This file is used to create count file for do the differential gene expression analysis of 4 patient samples (tumor/normal)

RAW files should be downloaded from McGill epigenome mapping centre. and need to create a table to run the job using file names of samples

table name: table_rnaseq4pat.txt

inputs: all the BAM files from patient RL354, RL371, RL380 and RL400 (mapped in to GRCh19)
we need to remap them to GRCh38
in the server RL400 has fastq files (so for that files we need to make a seperate script

output: outputs are stored in RNA-seq_counts_4patients folder. It consist of count files of RNA-seq reads for each sample (both normal and tumor)
It can be used for DESeq2 analysis

##########################################################################################
##########################################################################################
##########################################################################################

Description about   RNA-seq_4patients_diff.Rmd

Differential gene expression analysis of 4 patients  (tumor /normal)

inputs: count files from previous file stored in RNA-seq_counts_4patients (GRCh 38)

output: ordered.down.regulated.genes.csv (GRCh 38)

File contains differential gene expression information of genes

Columns

"baseMean expression values","log2FoldChange","log fold change Standard Error","stat","pvalue","padj"

##########################################################################################
##########################################################################################
##########################################################################################

Description about wgbs_step1.Rmd

This file is used to create wgbs input files for methyl kit to do a differential methylation analysis and also it does the differentila methylation analysis

Input files were obtained from epignome mapping centre. they were processed using Bismark. They are standard bismark output coverage files (csv) (GRCh38)

we can prepare them for methyl kit using our script

then these files were used for the differential dna methylation analysis

output : all.diff.methylated.min.per.grp.3.cpgs.hg38.txt (GRCh38)

this file include differentially methylated information for CpGs in the genome
columns

seqnames =chromosome name
start = start position of the base pair   
end = end position of the base pair   
width (blank)   
strand  
pvalue = pvalue of the differential methylation analysis 
qvalue  = multiple hypothesis testing 
meth.diff = similar to beta value 

##########################################################################################
##########################################################################################
##########################################################################################

Description about job_H3K27ac.pbs

This file is used to create peak calling files of histone marks from fastq files obtained from epigenome mapping centre 

Relevent table files for file name is included but you have to modify the scrupt for other histone marks
(GRCh38)

MACS2 will be used to call peaks and those outputs will be used for other analysis























