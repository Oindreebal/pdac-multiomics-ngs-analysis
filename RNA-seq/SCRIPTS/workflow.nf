nextflow.enable.dsl=2

// ================= PARAMETERS =================
params.samplesheet = "samples.csv"
params.index  = "/run/media/oindree-bal/Dexton SSD/WORKDIR/RNASEQ/index/GRCh38"
params.gtf    = "/run/media/oindree-bal/Dexton SSD/WORKDIR/RNASEQ/gencode.v46.annotation.gtf"
params.output = "/run/media/oindree-bal/Dexton SSD/WORKDIR/RNASEQ/nextflow_output"

// Keep memory usage low
params.threads = 2

// ================= HISAT2 =================
process HISAT2 {

    tag "$sample_id"

    cpus 2
    memory '6 GB'
    maxForks 1

    publishDir "${params.output}/hisat2", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          path("${sample_id}.sam"),
          path("${sample_id}.hisat2_mapstats.txt")

    script:
    """
    hisat2 \
        -p ${task.cpus} \
        --very-sensitive \
        --rna-strandness RF \
        -x ${params.index} \
        -1 ${read1} \
        -2 ${read2} \
        -S ${sample_id}.sam \
        2> ${sample_id}.hisat2_mapstats.txt
    """
}

// ================= SORT BAM =================
process SamtoolsSort {

    tag "$sample_id"

    cpus 2
    memory '4 GB'
    maxForks 1

    publishDir "${params.output}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file), path(mapstats)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam_file} | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam

    rm -f ${sam_file}
    """
}

// ================= FEATURECOUNTS =================
process FeatureCounts {

    cpus 2
    memory '4 GB'
    maxForks 1

    publishDir "${params.output}/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "featureCounts.txt"
    path "featureCounts.txt.summary"

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${params.gtf} \
        -o featureCounts.txt \
        -p \
        --countReadPairs \
        -B \
        -C \
        -t exon \
        -g gene_id \
        ${bam_files.join(' ')}
    """
}

// ================= WORKFLOW =================
workflow {

    samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }

    aligned = HISAT2(samples)

    sorted = SamtoolsSort(aligned)

    bam_files = sorted.map { sample, bam, bai -> bam }.collect()

    FeatureCounts(bam_files)
}
