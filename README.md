# Statistical and Molecular Methods for the Attribution of Body Fluids to Individual Donors
This repository contains details for the data processing and analysis of RNA/cDNA and gDNA (genomic DNA) sequencing data for body fluid identification and donor attribution of body fluids to donors using cSNPs (coding region SNPs). 
All code examples are written in the context of processing RNA/cDNA FASTQ files, but the same code was also applied to gDNA data.
Created as part of a Master of Science (MSc) thesis in Forensic Science at the University of Auckland.
## Contents:
* [Quality checking](#quality-checking)
  * [FastQC](#quality-check-of-raw-sequences-with-fastqc)
  * [MultiQC](#combine-fastqc-analyses-into-single-report-with-multiqc)
  * [Adapter removal](#adapter-removal)
* [Mapping and counting reads](#mapping-and-counting-reads)
  * [Mapping reads with BWA](#map-reads)
  * [Processing SAM files](#process-sam-files)
  * [Get mapping statistics](#mapping-stats)
  * [Count alignments with featureCounts](#count-alignments)
* [Plot alignment data](#plot-alignment-data)
* [Principal Component Analysis (PCA)](#pca)
* [Variant (cSNP) calling](#variant-calling)
* [Genotype analysis](#genotype-analysis)
* [Bayesian Network analysis](#bayesian-network-analysis)
## Quality checking
### Quality check of raw sequences with FastQC
All FASTQ files (forward/R1 and reverse/R2) were transferred to a raw directory and processed with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
cd cDNA_raw
module purge
module load fastqc
fastqc -o QC/cDNA_Raw_1/*
```
### Combine FastQC analyses into single report with MultiQC
Aggregate results into a single report with [MultiQC](https://multiqc.info/)
```bash
module purge
module load multiqc
cd ~/cDNA_seq
mkdir multiqc
#Copy QC files into multiQC directory
cp QC/*.fastqc multiqc/ 
cd multiqc/
multiqc .
```
### Adapter removal
Clean reads with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
```bash
cd ~/cDNA_seq
mkdir trimmed
module purge
module load cutadapt
cd cDNA_raw
for filename in *.fastq
do 
base=$(basename ${filename} .fastq)
cutadapt -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ../trimmed/${base}.trimmed.fastq ${filename} > ../trimmed/${base}.log
done
cd ../Trimmed
#Copy cutadatpt log files into MultiQC folder and add cutadapt information to report
cd ../multiqc
cp ../Trimmed/*log .
multiqc .
```
## Mapping and counting reads
### Map reads
Align sequences to reference genome (GRCh38/hg38; available from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)) using [Burrow-Wheelers Alignment (BWA)](https://github.com/lh3/bwa)
```bash
cd ~/cDNA_seq
module purge
module load bwa
mkdir refgenome
mkdir align
#Upload, unzip and index reference genome
cd refgenome
gunzip hg38.fa.gz
bwa index hg38.fa
cd ~/cDNA_seq/trimmed
#Align to reference genome
for filename in *R1_001.trimmed.fastq
do
base=$(basename ${filename} R1_001.trimmed.fastq)
bwa mem ../refgenome/hg38.fa $filename ${base}R2_001.trimmed.fastq > ../align/${filename%_R1_001.trimmed.fastq}.paired.sam
done
```
### Process SAM files
Process SAM files with [SAMtools](https://www.htslib.org/)
```bash
cd ~/cDNA_seq/Align
module purge
module load samtools
#Convert SAM files to BAM files
for filename in *.sam
do
base=$(basename ${filename} .sam)
samtools view -S -b ${filename} -o ${base}.bam
done
#Mark duplicates and index
cd ~/cDNA_Seq
mkdir Samtools
cd Samtools/
cp ../align/*.paired.bam . 
for filename in *.paired.bam
do
base=$(basename ${filename} .paired.bam)
samtools collate -@ 4 -O -u ${filename} | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - ${base}.mrkdup.bam
samtools index ${base}.mrkdup.bam
done
```
### Mapping stats
Generate mapping statistics with [SAMtools](https://www.htslib.org/)
```bash
for filename in *.mrkdup.bam
do
base=$(basename ${filename} .mrkdup.bam)
samtools idxstats ${filename} | awk '$3>10' | sort -k3 -n -r > ${base}.idx_sort.txt
samtools stats ${filename} > ${base}.stats.txt
done
#Add to MultiQC report if desired
```
### Count alignments
Count read alignments with featureCounts ([Subread package](https://subread.sourceforge.net/)). Human GTF file is available from [Gencode](https://www.gencodegenes.org/human/)
```bash
#Download human GTF file into refgenome directory
cd ~/cDNA_seq/Samtools
mkdir subread
cd subread
#Download Subread software and install using instructions provided (https://subread.sourceforge.net/)
cp ../*.mrkdup.bam . 
subread-2.0.6-source/bin/featureCounts -a refgenome/gencode.v45.chr_patch_hapl_scaff.annotation.gtf -o ./cDNA_counts.txt -T 2 -t exon -g gene_id ./*mrkdup.bam -p
```
Process featureCounts alignment count file in R
```R
#Extract only rows with alignments
hits <- cDNA_counts[!apply(cDNA_counts[, 7:47]==0, 1, all),] 
write.csv(hits, "fc_cDNA", row.names = FALSE) #Create csv file
#Further processing in MS Excel
```
## Plot alignment data
Plot alignment data for each gene in R using [ggplot2](https://ggplot2.tidyverse.org/), and arrange plots with [ggpubr](https://rpkgs.datanovia.com/ggpubr/)
```R
#Load packages
library(ggplot2)
library(ggpubr)
#Read in summarised specificity data for each gene
df <- read.delim("~/R/Specificity_data.txt")
#Assign colours for marker types
colours <- c(Blood = "#C74F8A", Saliva = "#0563C1", Vaginal = "#477942", Skin = "#FFB000")
#Plot specificity data for blood samples with ggplot2
Blood <- ggplot(df) + geom_col(aes(Gene, Blood, fill=Marker_Type),    
	show.legend = FALSE) + ylab("Reads") + ggtitle("Blood") + 
	theme(axis.text.x = element_text(angle = 90), 
	plot.title = element_text(hjust = 0.5), 
	axis.title = element_text(face = 'italic')) + 
	scale_fill_manual(values=colours) + 
	scale_x_discrete(limits=df$Gene) + 
	scale_y_continuous(breaks = scales::pretty_breaks(n = 8), 
	labels = function(x) format(x, scientific = TRUE)) 
Blood
#Repeat for each sample type - saliva, vaginal, and skin
#Arrange plots together with ggpubr
ggarrange(Blood, Saliva, Skin, Vaginal, nrow = 2, ncol = 2)
```
Plot alignment data for each sample in R using [ggplot2](https://ggplot2.tidyverse.org/), cleaned with [tidyr](https://tidyr.tidyverse.org/) and [gtools](https://github.com/cran/gtools)
```R
#Load packages
library(ggplot2)
library(tidyr)
library(gtools)
library(ggpubr)
#Read in alignment files for each sample type
BL <- read.delim("~/R/Counts_BL.txt") #e.g., blood sample data
#Assign colours
colours <- c(Blood = "#C74F8A", Saliva = "#0563C1", Vaginal = "#477942", Skin = "#FFB000" )
#Clean data using tidyr
BL <- gather(BL, 'key', 'value', -Sample)
#Order sample names using gtools
BL <- BL[gtools::mixedorder(BL$Sample), ]
#Plot with ggplot2
Blood <- ggplot(BL, aes(Sample,  value, fill = key)) +
	geom_bar(stat = "identity") + 
	ylab("Reads")+ xlab("Sample") + ggtitle("Blood") + 
	theme(axis.text.x = element_text(angle = 90), 
	plot.title = element_text(hjust = 0.5), 
	axis.title = element_text(face = 'italic')) + 
	scale_y_continuous(breaks = scales::pretty_breaks(n = 8), 
	labels = function(x) format(x, scientific = TRUE))
	scale_fill_manual(values=colours) +  
Blood
#Repeat for each sample type - saliva, skin, vaginal, and controls)
#Arrange sample plots with ggpubr – controls plotted seperately 
ggarrange(Blood, Saliva, Skin, Vaginal, nrow = 2, ncol = 2)
```
## PCA
Principal component analysis (PCA) in R using [ggplot2](https://ggplot2.tidyverse.org/) and [ggfortify](https://github.com/sinhrks/ggfortify)
```R
#Read in file with alignment data
#Load packages
library(ggfortify)
library(ggplot2)
library(tidyverse)
#PCA
cDNA.pca <- prcomp(cDNA_PCA[,c(2:13)], center = TRUE, scale = TRUE) 
print(cDNA.pca)
summary(cDNA.pca)
#PCA plot 
plot <- autoplot(cDNA.pca, data = cDNA_PCA,colour = 'Type') + scale_color_manual(values = c("#C74F8A","#0563C1","#FFB000","#477942"))
plot
#PCA plot with loadings 
plot.load <- autoplot(cDNA.pca, data = cDNA_PCA,colour = 'Type',loadings=TRUE, loadings.label=TRUE, loadings.label.size=3.4, loadings.colour='black', loadings.label.hjust = 1.2, loadings.label.vjust = 0) + scale_color_manual(values = c("#C74F8A","#0563C1","#FFB000","#477942"))
plot.load
#PCA plot with point labels
plot.text <- plot + geom_text(vjust=-1, label=cDNA_PCA$Sample)
plot.text
##Scree plot
#Calculate variance
explained_variance <- data.frame(PC=paste0("PC",1:12),                                
var_explained=(cDNA.pca$sdev)^2/sum((cDNA.pca$sdev)^2)*100)
#Create scree bar plot
scree.col <- ggplot(explained_variance, aes(x=fct_inorder(PC), y=var_explained, group=1))+ geom_col() + xlab("Principal Component") + ylab("Variance Explained (%)")
scree.col
```
## Variant calling
Variant calling with [BCFtools](https://github.com/samtools/bcftools) mpileup tool
```bash
module purge
module load bcftools 
cd ~/cDNA_seq/Samtools
#Copy BAM files into new variants directory
mkdir Variants
cp *mrkdup.bam Variants/
cd Variants/
#Write list of all BAM files
for file in *.mrkdup.bam 
do
echo $file >> bam_files.txt
done
#Create VCF with allele depth information
#Collate all samples into single vcf at this point
bcftools mpileup --max-depth 5000 -a FORMAT/AD,FORMAT/DP,INFO/AD -f ../refgenome/hg38.fa -b bam_files.txt | bcftools call -mv -Ob -o Variants/RNA_calls.raw.bcf
bcftools view RNA_calls.raw.bcf > RNA_calls.raw.vcf
```
Extract and filter cSNPs with [VCFtools](https://vcftools.github.io/index.html), [BCFtools](https://github.com/samtools/bcftools), and [SAMtools](https://www.htslib.org/)
```bash
module purge
module load vcftools
module load bcftools
Module load samtools
cd ~/cDNA_seq/Samtools/Variants
#Remove InDels, keep only SNPs
bcftools view --types snps cDNA_all_annot.vcf > cDNA_all.vcf
#Filter on allelic depth 
bcftools filter -S . -i 'FMT/AD[*:*]>7' -O v -o RNA_filt.vcf cDNA_all.vcf
#Extract SNPs based on position in genome 
bgzip RNA_filt.vcf
bcftools index RNA_filt.vcf.gz
bcftools view -R position.txt RNA_filt.vcf.gz > cDNA_SNPs.vcf
#Get allele frequency data
vcftools --vcf cDNA_SNPs.vcf --freq --out cDNA_SNPs
```
## Genotype analysis 
Plot allele depth data for a given participant in R with [ggplot2](https://ggplot2.tidyverse.org/)
```R
#Load packages
library(ggplot2)
library(tidyr)
library(tidyverse)
#Read in file with allele depth information for each ALT and REF allele per sample sourced from a given participant
Allele_prop <- read.delim("~/R/P9_prop.txt")
#Set colours
colours <- c(REF = "#C895E4", ALT = "#679E80")
#Organise data with tidyr
P9 <- gather(Allele_prop, key = "measure", value = "value", c("Blood_A", "Blood_B", "Saliva_A", "Saliva_B", "Skin_A", "Skin_B", "Vaginal_A", "Vaginal_B"))
#Plot with ggplot2 
Plot <- ggplot(P9, aes(x=SNP.ID, y=value, fill = Position)) +
	geom_bar(stat = "identity", position = "fill") + 
	scale_fill_manual(values=colours) + 
	ylab("Allele depth (%)") + 
	xlab("cSNP ID") + ggtitle("Participant 9") + 
	theme(axis.text.x = element_text(angle = 90), 
	plot.title = element_text(hjust = 0.5), 
	axis.title = element_text(face = 'italic'), 
	text = element_text(size = 13)) + 
	scale_y_continuous(breaks = scales::pretty_breaks(n = 5), 
	labels = scales::percent) + facet_wrap(~measure, ncol = 2) + 
	guides(fill=guide_legend(title="Allele"))
Plot
```
Plot allele depth ratio versus total allele depth in R
```R
#Load packages
library(ggplot2)
#Read in file with reference allele ratios, total allele depth, and corresponding sample types and marker types for cDNA or gDNA
df <- read.delim("~/R/cDNA_ratio", header = TRUE)
#Assign colours
colours <- c(Blood = "#C74F8A", Saliva = "#0563C1", Vaginal = "#477942", Skin = "#FFB000")
#Plot with ggplot2
Plot <- ggplot(df, aes(x=Total, y=Ratio)) + 
	geom_point(aes(color = Sample_Type,shape=Marker_Type),
  size=2) + scale_color_manual(values=colours) + 
	ylab("Ratio reference allele") + 
	xlab("Total allele depth") + 
	theme(text = element_text(size = 15))
Plot
```

## Bayesian Network analysis
Plot LR results in R with [ggplot2](https://ggplot2.tidyverse.org/) and [scales](https://scales.r-lib.org/) for every sample and hypothesis combination
```R
#Load packages
library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
#Read in file with LR results for all blood samples
df <- read.delim("~/R/BL_LR")
#Organise data with tidyr
df2 <- gather(df, key = "measure", value = "value", c("Other", "Blood", "Saliva", "Skin", "Vaginal"))
#Remove NA values with tidyr
df2 <- drop_na(df2)
#Set order for categorical variables (blood, saliva, skin etc.)
df2 <- df2[with(df2,order(measure)),]
#Set colours
colours <- c(Blood = "#C74F8A", Saliva = "#0563C1", Vaginal = "#477942", Skin = "#FFB000", Other = "#C895E4")
#Plot using ggplot2, with horizontal line at LR = 1 and log10 Y axis
Plot_BL <- ggplot(df2) + geom_col(aes(Sample, value, fill=fct_inorder(measure)), position = "dodge") +
	ylab("Likelihood Ratio (LR)") + ggtitle("Hp: Blood") + theme(axis.text.x = element_text(angle = 90),
	plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + scale_y_log10(breaks = log_breaks()) +
	guides(fill=guide_legend(title="Hd")) + geom_hline(yintercept = 1) + scale_fill_manual(values=colours)
Plot_BL     
#Redo Y axis labels for log scale
BL2 <- Plot_BL + scale_y_log10(breaks = sort(c(as.numeric(ggplot_build(Blood)$layout$panel_params[[1]]$y.labels), 1, 10, 100, 1000, 10000, 100000, 1000000)))
BL2
#Repeat for each sample type (saliva, skin, and vaginal)
```
Plot verbal equivalents in R
```R
#Load packages
library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
#Read in file with verbal equivalent data for all blood samples
df <- read.delim("~/R/BL_verb")
#Organise data with tidyverse
df2 <- gather(df, key = "measure", value = "value", c("Other", "Saliva", "Skin", "Vaginal"))
#Set order for categorical variables (blood, saliva, skin etc.)
df2 <- df2[with(df2,order(measure)),]
#Set colours
colours <- c(Blood = "#C74F8A", Saliva = "#0563C1", Vaginal = "#477942", Skin = "#FFB000", Other = "#C895E4")
#Plot using ggplot2
BL_verb <- ggplot(df2) + geom_col(aes(fct_inorder(LR), value, fill=fct_inorder(measure)), position = "stack") +
	ylab("Number of LR values") +xlab("")+ ggtitle("Hp: Blood") + theme(plot.title = element_text(hjust = 0.5),
	text = element_text(size = 15)) + guides(fill=guide_legend(title="Hd")) + scale_fill_manual(values=colours) +
	coord_flip()
BL_verb
#Repeat for each sample type
```

