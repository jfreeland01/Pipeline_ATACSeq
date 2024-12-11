# Pipeline_ATACSeq

## ** UNDER CONSTRUCTION **

Contact = jackfreeland01@gmail.com 

### **Outline** 

- [Introduction](#introduction)
- [Quality Control on Raw Reads](#quality-control-on-raw-reads)
- [Adapter Trimming](#adapter-trimming)
- [Quality Control on Trimmed Reads](#quality-control-on-trimmed-reads)
- [Alignment](#alignment)

## **Introduction**
This repository provides an example pipeline for processing **bulk ATAC-seq** (Assay for Transposase-Accessible Chromatin sequencing) data. The pipeline generates a consensus peak count matrix from raw fastq read files, followed by differential peak accessibility, overall TSS accessibility, and MOTIF enrichment analyses. This walkthrough will cover each step in detail, including the rational, code, and explanation of parameters/arguments when appropriate. An excellent resource to get familiarized with ATAC-seq analysis can be found in [Yan et al. (2020)](#https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-020-1929-3.pdf).

## **Quality Control on Raw Reads**
Before initiating formal processing steps, it is advisable to assess the overall quality of the raw sequencing files. This pipeline uses [FASTQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) to provide information on sequence quality, GC content, length distribution, duplicate sequences, overrepresented sequences, K-mer content and adapter contamination. In ATAC-seq data, it is expected to see Nextera transposase sequencing adapters over-represented and for the overall sequence quality to decrease near the 3' end. For pair-end sequencing, this should be run on both files. If interested, alternative/similar software  to [FASTQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are [FastP](#) or [BBDuk](https://sourceforge.net/projects/bbmap/).

```
# -t number of cores to be used
# -o output directory

fastqc -o <output_dir> <sample_ID>_R1.fastq.gz -t <#_of_cores> 

fastqc -o <output_dir> <sample_ID>_R2.fastq.gz -t <#_of_cores> 
```

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png" alt="Figure 1: FASTQC Adapter Sequence" width="600"/>

**Figure 1: Example FASTQC Output of Adapter Content**

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_SeqQual.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 2: Example FASTQC Output of Sequence Quality**







## **Adapter Trimming**
## **Quality Control on Trimmed Reads**
## **Alignment**


