# Pipeline_ATACSeq

## ** UNDER CONSTRUCTION **


## **Outline** 
### Background
- [Contact](#contact)
- [Introduction](#introduction)
### Pre-processing
- [Quality Control on Raw Reads](#quality-control-on-raw-reads)
- [Adapter Trimming](#adapter-trimming)
- [Quality Control on Trimmed Reads](#quality-control-on-trimmed-reads)
- [Alignment](#alignment)
- [Post-Alignment Filtering](#post-alignment-filtering)

### Post-processing
- [Peak Calling](#peak-calling)
- [Peak Count Matrix](#peak-count-matrix)
- [Differential Peak Analyses](#differential-peak-analyses)
- [MOTIF Enrichment](#motif-enrichment)
- [File Conversion Wig/bigWig](#file-conversion-wigbigwig)
- [Overall TSS Accessibility](#overall-tss-accessibility)

## **Contact**
For questions, comments, suggestions, anything, feel free to contact via git or through the following.

- Email: jackfreeland01@gmail.com
- LinkedIn: [@JackFreeland](https://www.linkedin.com/in/jack-freeland-384526142)
- Twitter: [@JackFreelandLab](https://x.com/JackFreelandLab)


## **Introduction**
This repository provides an example pipeline for processing **bulk ATAC-seq** data (Assay for Transposase-Accessible Chromatin sequencing). Starting with raw fastq files, the pipeline  calls individual peaks and can generate consensus peak count matrices, followed by differential peak accessibility, overall TSS accessibility, and MOTIF enrichment analyses. This walkthrough will cover each step in detail, including the rationale, example code, and explanations of parameters/arguments when appropriate. An excellent resource to first get familiarized with ATAC-seq analysis can be found in [Yan et al. (2020)](https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-020-1929-3.pdf).

## **Quality Control on Raw Reads**
Before initiating any formal processing steps, it is advisable to assess the overall quality of the raw sequencing files. This pipeline uses [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) that provide information on sequence quality, GC content, length distribution, duplicate sequences, overrepresented sequences, K-mer content and adapter contamination. For paired-end sequencing, this analysis should be performed on both files. In ATAC-seq [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf), we expect to see Nextera transposase sequencing adapters over-represented and a decrease in overall sequence quality near the 3' end. If interested, alternative or similar software to [FastQC](#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) includes [FastP](https://github.com/OpenGene/fastp) and [BBDuk](https://sourceforge.net/projects/bbmap/).

```
# -t    number of threads
# -o    output directory

fastqc -o <output_dir> <sample_ID>_R1.fastq.gz -t <#_of_threads>

fastqc -o <output_dir> <sample_ID>_R2.fastq.gz -t <#_of_threads>
```

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter.png" alt="Figure 1: FASTQC Adapter Sequence" width="600"/>

**Figure 1: Example FastQC Output of Adapter Content**

<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_SeqQual.png" alt="Figure 2: FASTQC Sequence Quality" width="600"/>

**Figure 2: Example FastQC Output of Sequence Quality**

## **Adapter Trimming**
As identified by [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), your reads are often contaminated by adapter sequences, which need to be trimmed before alignment. These adapter sequences are artificial sequences introduced during library preparation. If not removed, they can falsely align to the reference genome, leading to incorrect mapping results and increased noise. This pipeline uses [CutAdapt](https://cutadapt.readthedocs.io/en/stable/). If interested, alternative or similar software to [CutAdapt](https://cutadapt.readthedocs.io/en/stable/) includes [FastP](https://github.com/OpenGene/fastp), [ApadpterRemoval](https://adapterremoval.readthedocs.io/en/2.3.x/) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

Common adapter sequences are the following:

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

After running [CutAdapt](https://cutadapt.readthedocs.io/en/stable/), it is good practice to confirm the adapter sequences are no longer present by running [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) again. 

```
# -t    number of threads
# -o    output directory

fastqc -o <output_dir> <sample_ID>_R1.trim.fastq.gz -t <#_of_threads>

fastqc -o <output_dir> <sample_ID>_R2.trim.fastq.gz -t <#_of_threads>
```
<img src="https://github.com/jfreeland01/Pipeline_ATACSeq/blob/main/Figures/FASTQC_Adapter_trim.png" alt="Figure 3: FASTQC Adapter Sequence Post Cutadapt" width="600"/>

**Figure 3: Example FastQC Output of Adapter Content Post CutAdapt**

## **Alignment**
After confirming that adapter sequences have been trimmed, the next step is t align the reads to a reference genome. This pipeline uses [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to the human genome (GRCh38). [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) references can be found [here](https://benlangmead.github.io/aws-indexes/bowtie). In this example, GRCh38_noalt_decoy_as is used. 

The GRCh38 reference genome includes alternate haplotypes, which can cause reads to map equally to multiple regions, potentially resulting in lower quality scores. To address this, it is advantageous to remove alternate haplotypes and retain only the primary assembly ('noalt'). Additionally, including decoy sequences in your analysis can enhance accuracy. If a read aligns to a decoy sequence better than anywhere in the primary reference, it prevents false positives from misalignments within the primary genome. Both steps also decrease processing overhead and runtime. If interested, alternative or similar software to [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is [BWA](https://bio-bwa.sourceforge.net).

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
After aligning the reads, the data can be prepared for downstream analyses by converting the file type and applying QC/filtering measures:

- [Converting from SAM to Bam](#convert-from-sam-to-bam)
- [Sort Reads and Index](#sort-reads-and-index)
- [Adjust for Transposase Binding Offset](#adjust-for-transposase-binding-offset)
- [Filter Unaligned Reads](#filter-unaligned-reads)
- [Filter ChrM, DNA scaffold, and 'Decoy' Reads](#filter-chrm-dna-scaffold-and-decoy-reads)
- [Filter PCR Duplicates](#filter-pcr-duplicates)
- [Filter for Quality Reads](#filter-for-quality-reads)
- [Filter of ENCODE blacklists](#filter-of-encode-blacklists)

### **Convert from SAM to BAM**
Before applying filtering or QC, the SAM file (a text based format) generated by [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is converted to a BAM file (a binary format) using [Samtools](https://github.com/samtools/samtools). This conversion significant reduces storage requirements and processing time, as most tools are optimized to handle BAM files.
```
samtools view -@ <#_of_threads> \
-bS <sample_ID>.sam > <sample_ID>.bam
```

### **Sort Reads and Index**
Downstream software requires the BAM file to be sorted and indexed. Sorting organizes the alignments by their genomic coordinates, while indexing creates a companion index file that enables efficient random access to specific regions of the genome within the BAM files. Both steps are performed using [Samtools](https://github.com/samtools/samtools).
```
# -@    number of threads
# -o    output bam file

samtools sort \
-@ <#_of_threads> \
-o <sample_ID>_V1.bam \
<sample_ID>.bam

samtools index <sample_ID>_V1.bam
```

### **Adjust for Transposase Binding Offset**
Reads should be adjusted for transposase binding offset. Tn5 transposase binds to DNA and inserts sequencing adapters at staggered positions. It cuts the DNA with a 9-bp overhang between the strands, meaning the positions reported in sequencing data are not the exact sites of chromatin accessibility but are slightly shifted versions. Without adjustment, peak edges appear shifted, leading to imprecise peak summits. This offset can be corrected using [deepTools alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html).

```
# --bam     input BAM file
# -o        output BAM file

alignmentSieve --bam <sample_ID>_V1.bam \
-o <sample_ID>_V2.bam --ATACshift \
--numberOfProcessors <#_of_threads>
```

### **Filter Unaligned Reads**
Reads that failed to align to the reference genome should be removed as they do not provide any information about chromatin accessibility. Removing these reads also decreases file size and computational overhead. [Samtools](https://github.com/samtools/samtools) can be used to remove reads with the SAM flag 4, which indicates that the read is unmapped.

```
# -b    output results in BAM format
# -F    exclude reads with given flag

samtools view -@ <#_of_threads> -b -F 4 \ 
-o <sample_ID>_V3.bam <sample_ID>_V2.bam
```

### **Filter ChrM, DNA scaffold, and 'Decoy' Reads**
Mitochondiral reads should be removed as they, like unaligned reads, do not provide any information about chromatin accessibility (chrM). Unplaced and random DNA scaffolds should also be removed, as they are unlikely to reflect the chromatin accessibility of known, structured genomic loci. Decoy reads should be removed as decoy sequences are added to a reference genome to mitigate ambiguous read alignment. Decoy sequences largely do not correspond to function genomic loci and often represent technical artifacts or poorly chracterized regions. [Samtools](https://github.com/samtools/samtools) is used again.

```
# view      extracts BAM file and outputs as SAM (human-readable) to allow us to filter
# -@        number of threads
# -h        includes header (important for downstream)
# egrep -v  excludes lines containing the provided patern

# samtools view -b -o   conerts SAM back to BAM and writes file

samtools view -@ <#_of_threads> -h <sample_ID>_V3.bam | \
egrep -v "chrM|Un|random|decoy" | \
samtools view -b -o <sample_ID>_V4.bam
```

### **Filter PCR Duplicates**
PCR duplicates occur during library preparation. As the limited amount of DNA you start with in ATAC-seq requires PCR amplification to generate sufficient material for sequencing, multiple copies of the same DNA fragment are produced. This leads to identical (duplicated)reads that should be removed. To achieve this, the BAM file is first sorted again with [Samtools](https://github.com/samtools/samtools) and then filtered using [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard). 
```
samtools sort \
-@ <#_of_threads> \
-o <sample_ID>_V5.bam <sample_ID>_V4.bam

# -Xmx3G                  sets max memory allocation for the java virtual machine (JVM) 
# -XX:ParallelGCThreads   sets number of threads for java's garbage collection process
# -jar picard.jar         run picard.jar which contains Picard tools
# -M                      writes duplication metrics to file
# -I                      input file
# -O                      output file
# -REMOVE_DUPLICATES      if set to TRUE, duplicates are removed instead of just flagged

java -Xmx3G -XX:ParallelGCThreads=<#_of_threads> \
-jar picard.jar MarkDuplicates \
-M <sample_ID>_m.txt" \
-I <sample_ID>_V5.bam \
-O <sample_ID>_V6.bam \
-REMOVE_DUPLICATES true
```

### **Filter for Quality Reads**
Low-quality reads are then filtered out using [Samtools](https://github.com/samtools/samtools) to ensure downstream analyses are based on reliable and biologically meaningful data. 
```
# -q    MAPQ score threshold (higher = more strict)

samtools view -b -q 30 -@ <#_of_threads> \
<sample_ID>_V6.bam > <sample_ID>_V7.bam
```

### **Filter of ENCODE blacklists**
Finaly, the [ENCODE blacklist](https://www.nature.com/articles/s41598-019-45839-z) contains a list of genomic regions that are critical to remove when analyzing functional genomic data. Reads that map to these regions are typically not due to true biological signal but rather technical artifacts (e.g., misalignment, PCR amplification biases). Keeping these reads can introduce noise and false positives. [Bedtools](https://bedtools.readthedocs.io/en/latest/) can be used to filter out these regions. The final BAM file will also be indexed to aid in downstream analyses.
```
# subtract     subcommand used to subtract regions from one file based on overlaps with another
# -a    BAM file to be filtered
# -b    ENCODE blacklist (BED file)
# -A    only remove reads if they completely overlap

bedtools subtract -a <sample_ID>_V7.bam -b <hg38-blacklist.bed> \
-A > <sample_ID>_V8.bam
```

## **Peak Calling** ##
After filtering the BAM files, peaks can be called. Peaks are regions of the genome with high levels of read enrichment, indicating areas of open chromatin or accessible DNA. These regions are typically associated with regulatory elements, such as promoters, enhancers, transcription factor binding sites, or other DNA elements where the chromatin is less compact, allowing transcriptional machinery and regulatory proteins to bind. This pipeline uses [MACS3 callpeak](https://github.com/macs3-project/MACS) to identify peaks.

```
# -f    Specify format of input file (BAMPE = paired-end BAM files)
# -g    Specify gemoe size (hs = human genome, mm = mouse)
# -q    False discovery rate (FDR) threshold for peak calling
# -t    Path to input BAM file
# -n    Specifies a name prefix for output files

macs3 callpeak \
-f BAMPE \
-g hs \
-q 0.01 \
-t <sample_ID>_V8.bam \
-n <sample_ID> \
--outdir <output_directory>
```

## **Peak Count Matrix** ##
After calling peaks, a peak count matrix can be generated for downstream analyses. As a peak count matrix is a tabular representation of read counts for identified genomic regions (peaks) across multiple samples, it is essential to first generate a set of consensus peaks across all samples. This is because any specific region or peak representing the same biological element across samples may vary by a few bases due to differences in MACS3 peak calling or biology. Such variations can complicate downstream quantification, analysis, and interpretation of results.

### **Generate Concensus Peak File**
```
#!/bin/bash

### Find all narrowPeak files from MACS3 and concatenate them into one file
# $peak_dir     Directory containing all MACS3 output files

find "$peak_dir" -name "*.narrowPeak" -exec cat {} + > "$peak_dir/all_concatenate.narrowPeak.bed"


### Sort the combined file

sort -k1,1 -k2,2n "$peak_dir/all_concatenate.narrowPeak.bed" > "$peak_dir/all_concatenate_sorted.narrowPeak.bed"


### Merge the peaks

bedtools merge -i "$peak_dir/all_concatenate_sorted.narrowPeak.bed" > "$peak_dir/all_concatenate_sorted_merged.narrowPeak.bed"
```
### **Generate Concensus Count File**
```
#!/bin/bash

###  Create an array of all BAM files
# $bam_dir     Directory containing all processed BAM files
bam_files=($bam_dir/*V8.bam)

### Create header with BAM filenames for the multicov output, ensuring tab separation

header="chr\tstart\tend"
for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" | sed 's/_V8.bam//g')
    header="$header\t$sample_name"
done

### Output the header to the file
echo -e "$header" > "$peak_dir/multicov_merged_counts.txt"

### Run bedtools multicov to count reads in merged peaks across all BAM files
bedtools multicov -bams "${bam_files[@]}" -bed "$peak_dir/all_concatenate_sorted_merged.narrowPeak.bed" >> "$peak_dir/multicov_merged_counts.txt"
```

## **Differential Peak Analyses**
To perform differential peak analyses, workflows very similar to those used for differential gene expression analyses can be applied, including packages such as [DESeq2](#https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](#https://bioconductor.org/packages/release/bioc/html/edgeR.html), and [limma](#https://bioconductor.org/packages/release/bioc/html/limma.html).

Of note, the 'Peak Count Matrix' protocol will generate a count file with the chromosome number, starting base pair position, and ending base pair position in three separate columns. These columns will need to be combined into a single 'name' column before running any of the above packages. For example, by running DESeq2, you will generate a table as follows:

```
                    baseMean    log2FC      lfcSE       stat        pvalue      padj
chr1_808161_808284  30.10598    -0.380062   0.417783    -0.909712   0.362974    0.999999
chr1_817195_817520  84.26785    0.3425648   0.254628    1.3453500   0.178512    0.999999
...                 ...         ...         ...         ...         ...         ...
```
## **MOTIF Enrichment**
After performing differential peak analysis, motif enrichment analysis allows you to identify transcription factor binding sites or regulatory elements that are enriched within the differentially accessible regions. This provides insights into the potential regulatory mechanisms driving changes in chromatin accessibility and gene expression, helping to connect observed epigenetic changes with underlying biological processes or pathways. By uncovering enriched motifs, you can prioritize key transcription factors or regulatory networks for further functional validation. This pipeline uses [HOMER](#http://homer.ucsd.edu/homer/motif/) to perform the enrichment.  

## **File Conversion Wig/bigWig** ##
Many packages which visualize genomic data (such as ATAC) requires BAM files to be converted to either WIG or BigWig files. Here, [deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) is used to convert from BAM to BigWig and [UCSC Genome Broswer bigWigtoWig](https://www.encodeproject.org/software/bigwigtowig/) to convert from BigWig to Wig.

```
findMotifsGenome.pl <significant_peaks_file> hg38 <output_dir> -size given -p <#_of_CPUs>
```

### **BAM to BigWig**
```
# -b                        Input BAM file
# --normalizeUsing          Normalization methods to account for sequencing depth, genome size, or scaling factors.
# --effectiveGenomeSize     Specify effective genome size (portion of the genome mappable by reads)
# -p                        Number of threads for parallel processing
# -o                        Output file path
# --binSize                 Size of the bins (in base pairs) used to compute coverage

bamCoverage -b <sample_ID>_V8.bam \
--normalizeUsing RPGC \
--effectiveGenomeSize 2913022398 \ (for homo sapiens)
-p <#_of_threads> \
-o <sample_ID>.bw \
--binSize 10
```
### **BigWig to WIG**

```
bigWigToWig <sample_ID>.bw <GRCH38_noalt_decoy_as.chrom.sizes> <sample_ID>.wig

# In this example, the reference file <GRCH38_noalt_decoy_as.chrom.sizes> is structured as follows:

# chr1    248956422
# chr2    242193529
# chr3    198295559
# ...     ...
```

## **Overall TSS Accessibility** ##
