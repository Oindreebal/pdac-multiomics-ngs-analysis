#!/bin/bash

set -e

############################
# CONFIGURATION
############################

THREADS=2

WORKDIR="/home/oindree-bal/RNASEQ"
INDEX="$WORKDIR/index/GRCh38_index"
GTF="$WORKDIR/gencode.v46.annotation.gtf"
OUTDIR="$WORKDIR/results"

mkdir -p "$OUTDIR"/{hisat2,bam,counts}

#!/bin/bash

set -euo pipefail

#########################################
# DIRECTORIES
#########################################

THREADS=8

WORKDIR="/path/to/your/workdir"
OUTDIR="$WORKDIR/outputs"

mkdir -p "$OUTDIR/fastqc/raw"
mkdir -p "$OUTDIR/fastqc/trimmed"
mkdir -p "$OUTDIR/trimmed"
mkdir -p "$OUTDIR/logs"

# Trimmomatic adapter file
ADAPTERS="/path/to/Trimmomatic/adapters/TruSeq3-PE.fa"


#########################################
# QC NORMAL SAMPLE
#########################################

echo "========================================="
echo " FASTQC - NORMAL RAW READS"
echo "========================================="

fastqc \
    -t "$THREADS" \
    -o "$OUTDIR/fastqc/raw" \
    "$WORKDIR/Samples/Normal_1.fastq.gz" \
    "$WORKDIR/Samples/Normal_2.fastq.gz"


echo "========================================="
echo " TRIMMOMATIC - NORMAL"
echo "========================================="

trimmomatic PE \
    -threads "$THREADS" \
    -phred33 \
    "$WORKDIR/Samples/Normal_1.fastq.gz" \
    "$WORKDIR/Samples/Normal_2.fastq.gz" \
    "$OUTDIR/trimmed/Normal_1_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Normal_1_unpaired.fq.gz" \
    "$OUTDIR/trimmed/Normal_2_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Normal_2_unpaired.fq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    > "$OUTDIR/logs/Normal.trimmomatic.log" 2>&1


echo "========================================="
echo " FASTQC - NORMAL TRIMMED READS"
echo "========================================="

fastqc \
    -t "$THREADS" \
    -o "$OUTDIR/fastqc/trimmed" \
    "$OUTDIR/trimmed/Normal_1_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Normal_2_trimmed.fq.gz"


#########################################
# QC TUMOR SAMPLE
#########################################

echo "========================================="
echo " FASTQC - TUMOR RAW READS"
echo "========================================="

fastqc \
    -t "$THREADS" \
    -o "$OUTDIR/fastqc/raw" \
    "$WORKDIR/Samples/Tumor_1.fastq.gz" \
    "$WORKDIR/Samples/Tumor_2.fastq.gz"


echo "========================================="
echo " TRIMMOMATIC - TUMOR"
echo "========================================="

trimmomatic PE \
    -threads "$THREADS" \
    -phred33 \
    "$WORKDIR/Samples/Tumor_1.fastq.gz" \
    "$WORKDIR/Samples/Tumor_2.fastq.gz" \
    "$OUTDIR/trimmed/Tumor_1_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Tumor_1_unpaired.fq.gz" \
    "$OUTDIR/trimmed/Tumor_2_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Tumor_2_unpaired.fq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    > "$OUTDIR/logs/Tumor.trimmomatic.log" 2>&1


echo "========================================="
echo " FASTQC - TUMOR TRIMMED READS"
echo "========================================="

fastqc \
    -t "$THREADS" \
    -o "$OUTDIR/fastqc/trimmed" \
    "$OUTDIR/trimmed/Tumor_1_trimmed.fq.gz" \
    "$OUTDIR/trimmed/Tumor_2_trimmed.fq.gz"


echo "========================================="
echo " PRE-ALIGNMENT QC COMPLETE"
echo "========================================="

#########################################
# ALIGNING NORMAL SAMPLE
#########################################

echo "===== Aligning NORMAL ====="

hisat2 \
-p $THREADS \
--very-sensitive \
--rna-strandness RF \
-x "$INDEX" \
-1 "$WORKDIR/Samples/trimmed/Normal_1_trimmed.fq.gz" \
-2 "$WORKDIR/Samples/trimmed/Normal_2_trimmed.fq.gz" \
2> "$OUTDIR/hisat2/Normal.hisat2.log" \
| samtools sort \
-@ $THREADS \
-m 512M \
-O BAM \
-o "$OUTDIR/bam/Normal.sorted.bam" -

samtools index "$OUTDIR/bam/Normal.sorted.bam"

#########################################
# ALIGNING TUMOR SAMPLE
#########################################

echo "===== Aligning TUMOR ====="

hisat2 \
-p $THREADS \
--very-sensitive \
--rna-strandness RF \
-x "$INDEX" \
-1 "$WORKDIR/Samples/trimmed/Tumor_1_trimmed.fq.gz" \
-2 "$WORKDIR/Samples/trimmed/Tumor_2_trimmed.fq.gz" \
2> "$OUTDIR/hisat2/Tumor.hisat2.log" \
| samtools sort \
-@ $THREADS \
-m 512M \
-O BAM \
-o "$OUTDIR/bam/Tumor.sorted.bam" -

samtools index "$OUTDIR/bam/Tumor.sorted.bam"

#########################################
# FEATURECOUNTS
#########################################

echo "===== Running FeatureCounts ====="

featureCounts \
-T $THREADS \
-a "$GTF" \
-p \
--countReadPairs \
-B \
-C \
-t exon \
-g gene_id \
-o "$OUTDIR/counts/featureCounts.txt" \
"$OUTDIR/bam/Normal.sorted.bam" \
"$OUTDIR/bam/Tumor.sorted.bam"

echo
echo "RNA-seq pipeline completed successfully."
