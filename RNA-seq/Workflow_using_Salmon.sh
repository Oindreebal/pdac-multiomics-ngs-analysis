#!/bin/bash

set -e

#########################################
# CREATE OUTPUT DIRECTORIES
#########################################

mkdir -p "$WORKDIR/results"/{fastqc_raw,fastqc_trimmed}
mkdir -p "$WORKDIR/Samples/trimmed"

#########################################
# CONFIGURATION
#########################################

THREADS=4

WORKDIR="/run/media/oindree-bal/Dexton SSD/RNASEQ"
FASTQ_DIR="$WORKDIR/Samples/trimmed"

TRANSCRIPTOME="$WORKDIR/ref/gencode.v46.transcripts.fa.gz"

GTF="$WORKDIR/ref/gencode.v46.annotation.gtf"

INDEX="$WORKDIR/index/salmon_index"

OUTDIR="$WORKDIR/results/salmon"

mkdir -p "$INDEX"
mkdir -p "$OUTDIR"


#########################################
# FASTQC BEFORE TRIMMING
#########################################

echo "===== FastQC : Raw Reads ====="

fastqc \
-t $THREADS \
-o "$WORKDIR/results/fastqc_raw" \
"$WORKDIR/Samples/Normal_1.fastq.gz" \
"$WORKDIR/Samples/Normal_2.fastq.gz" \
"$WORKDIR/Samples/Tumor_1.fastq.gz" \
"$WORKDIR/Samples/Tumor_2.fastq.gz"

#########################################
# TRIMMOMATIC
#########################################

TRIMMOMATIC=/usr/share/java/trimmomatic.jar

ADAPTERS=/usr/share/trimmomatic/TruSeq3-PE.fa

echo "===== Trimming NORMAL ====="

java -jar $TRIMMOMATIC PE \
-threads $THREADS \
-phred33 \
"$WORKDIR/Samples/Normal_1.fastq.gz" \
"$WORKDIR/Samples/Normal_2.fastq.gz" \
"$WORKDIR/Samples/trimmed/Normal_1_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Normal_1_trimmed_unpaired.fq.gz" \
"$WORKDIR/Samples/trimmed/Normal_2_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Normal_2_trimmed_unpaired.fq.gz" \
ILLUMINACLIP:${ADAPTERS}:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:36

echo "===== Trimming TUMOR ====="

java -jar $TRIMMOMATIC PE \
-threads $THREADS \
-phred33 \
"$WORKDIR/Samples/Tumor_1.fastq.gz" \
"$WORKDIR/Samples/Tumor_2.fastq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_1_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_1_trimmed_unpaired.fq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_2_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_2_trimmed_unpaired.fq.gz" \
ILLUMINACLIP:${ADAPTERS}:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:36

#########################################
# FASTQC AFTER TRIMMING
#########################################

echo "===== FastQC : Trimmed Reads ====="

fastqc \
-t $THREADS \
-o "$WORKDIR/results/fastqc_trimmed" \
"$WORKDIR/Samples/trimmed/Normal_1_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Normal_2_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_1_trimmed_paired.fq.gz" \
"$WORKDIR/Samples/trimmed/Tumor_2_trimmed_paired.fq.gz"

echo "PRE-PROCESSING COMPLETE"

#########################################
# BUILD SALMON INDEX (only first time)
#########################################

if [ ! -d "$INDEX" ] || [ -z "$(ls -A "$INDEX")" ]; then

    echo "Building Salmon index..."

    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX" \
        -k 31 \
        --threads $THREADS

fi

#########################################
# NORMAL
#########################################

echo "Quantifying NORMAL..."

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$FASTQ_DIR/Normal_1_trimmed.fq.gz" \
    -2 "$FASTQ_DIR/Normal_2_trimmed.fq.gz" \
    --validateMappings \
    -p $THREADS \
    -o "$OUTDIR/Normal"

#########################################
# TUMOR
#########################################

echo "Quantifying TUMOR..."

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$FASTQ_DIR/Tumor_1_trimmed.fq.gz" \
    -2 "$FASTQ_DIR/Tumor_2_trimmed.fq.gz" \
    --validateMappings \
    -p $THREADS \
    -o "$OUTDIR/Tumor"

#########################################
# CREATE tx2gene
#########################################

echo "Creating tx2gene..."

awk '
BEGIN{OFS="\t"}
/transcript_id/ && /gene_id/{
    match($0,/transcript_id "([^"]+)"/,a)
    match($0,/gene_id "([^"]+)"/,b)
    if(a[1]!="" && b[1]!="")
        print a[1],b[1]
}
' "$GTF" | sort -u > "$OUTDIR/tx2gene.txt"

#########################################
# tximport
#########################################

cat > "$OUTDIR/tximport.R" <<EOF
library(tximport)

samples <- c("Normal","Tumor")

files <- file.path(
    "$OUTDIR",
    samples,
    "quant.sf"
)

names(files) <- samples

tx2gene <- read.table(
    "$OUTDIR/tx2gene.txt",
    header=FALSE,
    sep="\t"
)

colnames(tx2gene) <- c("TXNAME","GENEID")

txi <- tximport(
    files,
    type="salmon",
    tx2gene=tx2gene
)

write.csv(txi\$counts,
          "$OUTDIR/gene_counts.csv")

write.csv(txi\$abundance,
          "$OUTDIR/gene_TPM.csv")
EOF

Rscript "$OUTDIR/tximport.R"

echo
echo "====================================="
echo "Pipeline completed successfully."
echo "====================================="
echo
echo "Outputs:"
echo "$OUTDIR/Normal"
echo "$OUTDIR/Tumor"
echo "$OUTDIR/gene_counts.csv"
echo "$OUTDIR/gene_TPM.csv"
