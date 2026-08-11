set -euo pipefail

###########################################
# SAMPLE ID
###########################################

TUMOR="ERR875340"

###########################################
# PATHS
###########################################
WORKDIR="/run/media/debojyoti/Dexton SSD/WORKDIR"
REF="$WORKDIR/ref/Homo_sapiens_assembly38.fasta"
DBSNP="$WORKDIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="$WORKDIR/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="$WORKDIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
GATK="$WORKDIR/gatk-4.5.0.0/gatk"
FUNCOTATOR="$WORKDIR/funcotator_dataSources.v1.8.hg38.20230908s"
OUTPUT="$WORKDIR/outputs"
RAW_VCF="$OUTPUT/vcf_files/${TUMOR}_somatic_raw.vcf.gz"
FILTERED_VCF="$OUTPUT/vcf_files/${TUMOR}_somatic_filtered.vcf.gz"
ANNOTATED_VCF="$OUTPUT/vcf_files/${TUMOR}_somatic_annotated.vcf.gz"

###########################################

mkdir -p \
"$OUTPUT/trimmed_reads" \
"$OUTPUT/bam_files" \
"$OUTPUT/logs" \
"$OUTPUT/qc_reports" \
"$OUTPUT/vcf_files"
    
###########################################
# TARGET PANEL
###########################################

TARGET_PANEL="$WORKDIR/targets.bed"

###########################################
# FastQC on raw reads
###########################################

echo "Running FastQC on raw reads..."

fastqc -t 8 \
    "$FASTQ_R1" \
    "$FASTQ_R2" \
    -o "$OUTPUT/qc_reports"

###########################################
# INPUT FASTQ
###########################################

FASTQ_R1="$WORKDIR/fastq/${TUMOR}_R1.fastq.gz"
FASTQ_R2="$WORKDIR/fastq/${TUMOR}_R2.fastq.gz"

mkdir -p \
"$OUTPUT/trimmed_reads" \
"$OUTPUT/bam_files" \
"$OUTPUT/logs"

###########################################
# Trimming
###########################################

echo "Trimming reads..."

trimmomatic PE -threads 8 \
    "$FASTQ_R1" \
    "$FASTQ_R2" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R1_unpaired.fq.gz" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R2_unpaired.fq.gz" \
    ILLUMINACLIP:"$WORKDIR/TruSeq3-PE.fa":2:30:10 \
    LEADING:3 TRAILING:3 MINLEN:36
    
###########################################
# FASTQC for trimmed reads
###########################################
    
fastqc -t 8 \
"$OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz" \
"$OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz" \
-o "$OUTPUT/qc_reports"

###########################################
# Alignment
###########################################

echo "Aligning reads..."

bwa mem -t 8 \
    -R "@RG\tID:${TUMOR}\tSM:${TUMOR}\tPL:ILLUMINA" \
    "$REF" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz" \
    "$OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz" \
| samtools sort -o "$OUTPUT/bam_files/${TUMOR}_sorted.bam"

###########################################
# Index BAM
###########################################

samtools index "$OUTPUT/bam_files/${TUMOR}_sorted.bam"

###########################################
# Mark Duplicates
###########################################

echo "Marking duplicates..."

"$GATK" MarkDuplicates \
    -I "$OUTPUT/bam_files/${TUMOR}_sorted.bam" \
    -O "$OUTPUT/bam_files/${TUMOR}_dedup.bam" \
    -M "$OUTPUT/logs/${TUMOR}_metrics.txt"

samtools index "$OUTPUT/bam_files/${TUMOR}_dedup.bam"

###########################################
# INPUT BAM
###########################################

INPUT_BAM="$OUTPUT/bam_files/${TUMOR}_dedup.bam"

RECAL_TABLE="$OUTPUT/logs/${TUMOR}_recal.table"
RECAL_BAM="$OUTPUT/bam_files/${TUMOR}_recal.bam"

###########################################
# BQSR (Targeted)
###########################################

echo "Running Base Quality Score Recalibration..."

"$GATK" BaseRecalibrator \
    -R "$REF" \
    -I "$INPUT_BAM" \
    --known-sites "$DBSNP" \
    --known-sites "$INDELS" \
    --known-sites "$MILLS" \
    -L "$TARGET_PANEL" \
    -O "$RECAL_TABLE"

echo "Applying BQSR..."

"$GATK" ApplyBQSR \
    -R "$REF" \
    -I "$INPUT_BAM" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -L "$TARGET_PANEL" \
    -O "$RECAL_BAM"

samtools index "$RECAL_BAM"

###########################################
# Somatic Variant Calling
###########################################

echo "Running Mutect2 on target panel..."

"$GATK" --java-options "-Xmx8g" Mutect2 \
    -R "$REF" \
    -I "$RECAL_BAM" \
    -tumor "$TUMOR" \
    -L "$TARGET_PANEL" \
    -O "$RAW_VCF"

###########################################
# Filter Variants
###########################################

echo "Filtering Mutect2 calls..."

"$GATK" FilterMutectCalls \
    -R "$REF" \
    -V "$RAW_VCF" \
    -O "$FILTERED_VCF"

###########################################
# Index VCF
###########################################

if [ ! -f "${FILTERED_VCF}.tbi" ]; then
    tabix -p vcf "$FILTERED_VCF"
fi

###########################################
# Funcotator Annotation
###########################################

echo "Annotating variants..."

"$GATK" Funcotator \
    --variant "$FILTERED_VCF" \
    --reference "$REF" \
    --ref-version hg38 \
    --data-sources-path "$FUNCOTATOR" \
    --output "$ANNOTATED_VCF" \
    --output-file-format VCF
echo "Target-panel variant calling completed."

