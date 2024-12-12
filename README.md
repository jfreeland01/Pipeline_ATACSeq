# Pipeline_ATACSeq

## ** UNDER CONSTRUCTION **

### **For questions, comments, suggestions, etc, contact**
- Email: jackfreeland01@gmail.com
- LinkedIn: [@JackFreeland](https://www.linkedin.com/in/jack-freeland-384526142)
- Twitter:  [@JackFreelandLab](https://x.com/JackFreelandLab)

### **Outline** 

- [Introduction](#introduction)
- [Quality Control on Raw Reads](#quality-control-on-raw-reads)
- [Adapter Trimming](#adapter-trimming)
- [Quality Control on Trimmed Reads](#quality-control-on-trimmed-reads)
- [Alignment](#alignment)
- [Post-Alignment Filtering](#post-alignment-filtering)

## **Introduction**
This repository provides an example pipeline for processing **bulk ATAC-seq** (Assay for Transposase-Accessible Chromatin sequencing) data. Starting with raw fastq reads, the pipeline  calls individual peaks and can generate consensus peak count matrices, followed by differential peak accessibility, overall TSS accessibility, and MOTIF enrichment analyses. This walkthrough will cover each step in detail, including the rational, example code, and explanations of parameters/arguments when appropriate. An excellent resource to first get familiarized with ATAC-seq analysis can be found in [Yan et al. (2020)](#https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-020-1929-3.pdf).

## **Quality Control on Raw Reads**
Before initiating formal processing steps, it is advisable to assess the overall quality of the raw sequencing files. This pipeline uses [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) to provide information on sequence quality, GC content, length distribution, duplicate sequences, overrepresented sequences, K-mer content and adapter contamination. For pair-end sequencing, this should be run on both files. In ATAC-seq [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf), we expect to see Nextera transposase sequencing adapters over-represented and for the overall sequence quality to decrease near the 3' end. If interested, alternative/similar software to [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are [FastP](#https://github.com/OpenGene/fastp) and [BBDuk](https://sourceforge.net/projects/bbmap/).

```
#   -t number of threads
#   -o output directory

fastqc -o <output_dir> <sample_ID>_R1.fastq.gz -t <#_of_threads>

fastqc -o <output_dir> <sample_ID>_R2.fastq.gz -t <#_of_threads>
```

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png" alt="Figure 1: FASTQC Adapter Sequence" width="600"/>

**Figure 1: Example FastQC Output of Adapter Content**

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_SeqQual.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 2: Example FastQC Output of Sequence Quality**

## **Adapter Trimming**
As was likely identified by [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), your reads are often contaminated by adapter sequences which need to be trimmed before alignment. These adapter sequences are artificial sequences introduced during library preparation. If not removed, these sequences can falsely align to the reference genome, leading to incorrect mapping results and increased noise. This pipeline uses [CutAdapt](#https://cutadapt.readthedocs.io/en/stable/). If interested, alternative/similar software to [CutAdapt](#https://cutadapt.readthedocs.io/en/stable/) are [FastP](#https://github.com/OpenGene/fastp), [ApadpterRemoval](#https://adapterremoval.readthedocs.io/en/2.3.x/) and [Trimmomatic](#http://www.usadellab.org/cms/?page=trimmomatic).

Common adapter sequences to be trimmed are the following:

- Illumina:   AGATCGGAAGAGC
- Small RNA:  TGGAATTCTCGG
- Nextera:    CTGTCTCTTATA
###

```
# -a    adapter sequence on the forward read (R1)
# -A    adapter sequence on the reverse read (R2)
# -j    number of threads
# -q    phred quality score threshold (q of 20 keeps bases with 99% accuracy)
# -O    minimum overlap between the adapter sequence and the read required for the adapter to be trimmed
# -m    mimimum read length rquired to be kept after trimming
# -o    output file for the trimmed forward reads (R1)
# -p    output file for the trimmed reverse reads (R2)

python3 -m cutadapt
-a CTGTCTCTTATA \ # Nextera example (most common)
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
#   -t number of threads
#   -o output directory

fastqc -o <output_dir> <sample_ID>_R1.fastq.gz -t <#_of_threads>

fastqc -o <output_dir> <sample_ID>_R2.fastq.gz -t <#_of_threads>
```
<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter_trim.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 3: Example FastQC Output of Adapter Content Post CutAdapt**

## **Alignment**
After confirming the adapter sequences were trimmed, we can now align the reads to a reference genome. This pipeline uses [Bowtie 2](#https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to the human genome (GRCh38). [Bowtie 2](#https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) references can be found [here](#https://benlangmead.github.io/aws-indexes/bowtie). In this example, GRCh38_noalt_decoy_as is used. 

Since the GRCh38 reference genome includes alternate haplotypes, it is advantageous to remove them and retain only the primary assembly ('noalt'). This prevents reads from mapping equally to multiple regions, which could result in lower quality scores. Including decoy sequences in your analysis can enhance accuracy. If a read aligns to a decoy sequence better than anywhere in the reference, anywhere that it would align within the reference without the decoys being included would result in a false positive. Both steps also decrease processing overhead and time. If interested, alternative/similar software to [Bowtie 2](#https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is [BWA](#https://bio-bwa.sourceforge.net).

```
# --non-deterministic   default, resolves ties randomly
# --mm                  enables memory mapping
# --phred33             default, FASTQ files follow Phred+33 scale
# --very-sensitive      slower but most sensitive alignment mode
# -p                    number of threads
# -x                    reference genome
# -1                    input file for the trimmed forward reads (R1)
# -2                    input file for the trimmed reverse reads (R2)
# -S                    output file in SAM format

bowtie2 --non-deterministic --mm --phred33 --very-sensitive \
-p <#_of_threads> \
-x <GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as> \
-1 <sample_ID>_R1.trim.fastq.gz \
-2 <sample_ID>_R2.trim.fastq.gz \
-S <sample_ID>.sam

```
## **Post-Alignment Filtering**
After aligning the reads, the data can now be prepped for downstream analyses by converting file type and aplying QC/filtering measures:

- [Converting from SAM to Bam](#convert-from-sam-to-bam)
- [Sort Reads and Index](#sort-reads-and-index)
- [Adjust for Transposase Binding Offset](#adjust-for-transposase-binding-offset)
- [Filter Unaligned Reads](#filter-unaligned-reads)
- [Filter ChrM and 'Decoy' Reads](#filter-chrm-and-decoy-reads)
- [Filter PCR Duplicates](#filter-pcr-duplicates)
- [Filter for Quality Reads](#filter-for-quality-reads)
- [Filter of ENCODE blacklists](#filter-of-encode-blacklists)

### **Convert from SAM to BAM**
Before filtering/QC, the SAM file (text based format) from [Bowtie 2](#https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is converted to a BAM file (binary format) using [Samtools](#https://github.com/samtools/samtools). This saves significant amounts of storage and processing time as most tools are optimized to process BAM files.
```
samtools view -@ <#_of_threads> \
-bS <sample_ID>.sam > <sample_ID>.bam
```

### **Sort Reads and Index**
Downstream software requires the BAM file to be sorted and indexed. Indexing creates a companion index file that allows for efficient random access to specific regions of the genome within the BAM files. This can all be done using [Samtools](#https://github.com/samtools/samtools).
```
# -@    number of threads
# -o    output bam file

samtools sort \
-@ <#_of_threads> \
-o <sample_ID>_V1.bam \
<sample_ID>.bam

samtools index <${sample_ID}_V1.bam>
```

### **Adjust for Transposase Binding Offset**
Reads should be adjusted for transposase binding offset. Tn5 transposase binds to DNA and inserts sequencing adapters at staggered positions. It cuts the DNA with a 9-bp overhang between the strands, meaning the positions it reports in sequencing data are not the precise positions of chromatin acessibility, but are slightly shifted version. Wihtout adjustment, peak edges appear shifted, leading to imprecise peak summits. This can be corrected using [deepTools alignmentSieve](#https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html).
```
alignmentSieve --bam <sample_ID>_V1.bam \
-o <sample_ID>_V2.bam --ATACshift \
--numberOfProcessors <#_of_threads>
```

### **Filter Unaligned Reads**
Reads that failed to to align to the reference genome should be removed as they do not provide any information about chromatin accessibility. Removing also decreases file size and computational overhead. [Samtools](#https://github.com/samtools/samtools) is used to remove reads with the SAM flag 4 indicating it is unmapped.

```
# -b    outputs results in BAM format
# -F    excludes reads with given flag

samtools view -@ <#_of_threads> -b -F 4 \ 
-o <sample_ID>_V3.bam <sample_ID>_V2.bam
```

### **Filter ChrM , scaffold, and 'Decoy' Reads**
Mitochondiral reads should be removed as they, like unaligned reads, do not provide any information about chromatin accessibility (chrM). Unplaced and random DNA scaffolds should be removed as they are unlikely to reflect chromatin accessibility of known, structured genomic loci. Decoy reads should be removed as decoy sequences are added to a reference genome to mitigate ambiguous read alignment. Decoy sequences largely do not correspond to function genomic loci and often represent technical artifacts or poorly chracterized regions. [Samtools](#https://github.com/samtools/samtools) is again used.


```
samtools view -@ <#_of_threads> -h <sample_ID>_V3.bam | \
egrep -v "chrM|Un|random|decoy" | \
samtools view -b -o <sample_ID>_V4.bam
```

### **Filter PCR Duplicates**
```
samtools sort -@ "$ncor" -o <sample_ID>_V5.bam <sample_ID>_V4.bam

java -Xmx3G -XX:ParallelGCThreads="$ncor" \
-jar "$ref_dir/picard.jar" MarkDuplicates \
-M "$bam_dir/dup_files/${sample_ID}_m.txt" \
-I <sample_ID>_V5.bam \
-O <sample_ID>_V6.bam \
-REMOVE_DUPLICATES true
```

### **Filter for Quality Reads**
```
samtools view -b -q 30 -@ <#_of_threads> <sample_ID>_V6.bam > <sample_ID>_V7.bam
```

### **Filter of ENCODE blacklists**
```
bedtools subtract -a <sample_ID>_V7.bam -b <hg38-blacklist.v2_sorted.bed> \
-A > <sample_ID>_V8.bam
```