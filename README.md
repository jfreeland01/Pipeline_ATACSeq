# Pipeline_ATACSeq

## *** UNDER CONSTRUCTION ***

Contact = jackfreeland01@gmail.com 

### Outline 

- [Introduction](#introduction)
- [Quality Control on Raw Reads](#quality-control-on-raw-read-files)
- [Adapter Trimming](#adapter-trimming)
- [Quality Control on Trimmed Reads](#quality-control-on-trimmed-reads)
- [Alignment](#alignment)

## Introduction

## Quality Control on Raw Read Files
Before initiating formal processing steps, it is advisable to assess the overall quality of the raw sequencing files. This pipeline uses [FASTQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) to provide information on sequence quality, GC content, length distribution, duplicate sequences, overrepresented sequences, K-mer content and adapter contamination. 

In ATAC-seq data expect to see Nextera transposase sequencing adapters over-represented and for the overall base sequence quality to fall off near the 3' end.

If interested, alternative/similiar softwares to [FASTQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are [FastP](#) or [BBDuk](https://sourceforge.net/projects/bbmap/).


<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png" alt="Figure 1: Adapter Sequence" width="500"/>


![Figure 1: FASTQC Adapter Sequence](https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png)



## Adaptor Trimming
## Quality Control on Trimmed Reads
## Alignment


