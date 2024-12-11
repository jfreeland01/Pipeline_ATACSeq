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
Before initiating formal processing steps, it is advisable to assess the overall quality of the raw sequencing files. This pipeline uses [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) to provide information on sequence quality, GC content, length distribution, duplicate sequences, overrepresented sequences, K-mer content and adapter contamination. In ATAC-seq data, it is expected to see Nextera transposase sequencing adapters over-represented and for the overall sequence quality to decrease near the 3' end. For pair-end sequencing, this should be run on both files. If interested, alternative/similar software  to [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are [FastP](#https://github.com/OpenGene/fastp) and [BBDuk](https://sourceforge.net/projects/bbmap/).

```
#   -t number of CPUs/threads to be used
#   -o output directory

fastqc -o <output_dir> <sample_ID>_R1.fastq.gz -t <#_of_CPUs> 

fastqc -o <output_dir> <sample_ID>_R2.fastq.gz -t <#_of_CPUs> 
```

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png" alt="Figure 1: FASTQC Adapter Sequence" width="600"/>

**Figure 1: Example FastQC Output of Adapter Content**

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_SeqQual.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 2: Example FastQC Output of Sequence Quality**

## **Adapter Trimming**
As identified by [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), your reads are often contaminated by adapter sequences which need to be trimmed before alignment. These adapter sequences are artificial sequences introduced during library preparation. If not removed, these sequences can falsely align to the reference genome, leading to incorrect mapping results and increased noise. This pipeline uses [CutAdapt](#https://cutadapt.readthedocs.io/en/stable/). If interested, alternative/similar software  to [CutAdapt](#https://cutadapt.readthedocs.io/en/stable/) are [FastP](#https://github.com/OpenGene/fastp), [ApadpterRemoval](#https://adapterremoval.readthedocs.io/en/2.3.x/) and [Trimmomatic](#http://www.usadellab.org/cms/?page=trimmomatic).

Common adapter sequences to be trimmed are the following:

- Illumina:   AGATCGGAAGAGC
- Small RNA:  TGGAATTCTCGG
- Nextera:    CTGTCTCTTATA
###

```
# -a    adapter sequence on the forward read (R1)
# -A    adapter sequence on the reverse read (R2)
# -j    number of CPUs/threads
# -q    sets a base phred quality score threshold (q of 20 keeps bases with 99% accuracy)
# -O    minimum overlap between the adapter sequence and the read required for the adapter to be trimmed
# -m    mimimum read length kept after trimming
# -o    output file for the trimmed forward reads (R1)
# -p    output file for the trimmed reverse reads (R2)

    python3 -m cutadapt
    -a CTGTCTCTTATA \ # Nextera example
    -A CTGTCTCTTATA \
    -j <#_of_CPUs> \
    -q 20 \
    -O 6 \
    -m 35 \
    -o <sample_ID>_R1.trim.fastq.gz \
    -p <sample_ID>_R2.trim.fastq.gz \
    <sample_ID>_R1.fastq.gz \
    <sample_ID>_R2.fastq.gz

```
## **Quality Control on Trimmed Reads**

After running [CutAdapt](#https://cutadapt.readthedocs.io/en/stable/), it is good practice to confirm the adapter sequences are no longer detected by running [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) again. 

```
#   -t number of CPUs/threads to be used
#   -o output directory

fastqc -o <output_dir> <sample_ID>_R1.trim.fastq.gz -t <#_of_CPUs> 

fastqc -o <output_dir> <sample_ID>_R2.trim.fastq.gz -t <#_of_CPUs> 
```
<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter_trim.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 3: Example FastQC Output of Adapter Content Post CutAdapt**

## **Alignment**


