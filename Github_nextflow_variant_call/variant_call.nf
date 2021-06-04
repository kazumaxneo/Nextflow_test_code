#!/usr/bin/env nextflow

/*
* Variant call jusing breseq
* Kazuma Uesaka
*/

PROGRAM_NAME = workflow.manifest.name

Channel
.fromFilePairs( params.reads )
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { read_pairs_ch }

// message
println "reads: $params.reads"


// process

/*
 * Step 1 cleaning
 */
process preprocessing {
    cpus 2
    tag "$pair_id"
    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    set pair_id, "${pair_id}_QT_R1.fq.gz" into read1_ch
    set pair_id, "${pair_id}_QT_R2.fq.gz" into read2_ch

    """
    fastp -i ${pair_id}_R1.fq.gz -I ${pair_id}_R2.fq.gz \
    -o ${pair_id}_QT_R1.fq.gz -O ${pair_id}_QT_R2.fq.gz \
    -h ${pair_id}_fastp.html -j ${pair_id}_fastp.json \
    -q 30 -n 5 -t 1 -T 1 -l 20 -w ${task.cpus}
    """
}

/*
 * Step 2 variant call
 */
process mapping {

    tag "$pair_id"
    publishDir params.outdir2, mode: 'copy'
    cpus 3
    input:
    path gff from params.gff
    path fasta from params.fasta
    tuple val(pair_id), path(reads1) from read1_ch
    tuple val(pair_id), path(reads2) from read2_ch

    output:
    set pair_id, "${pair_id}.bam" into bam_ch
    set pair_id, "${pair_id}*flagstats"

    """
    minimap2 -ax sr -t ${task.cpus} -R "@RG\\tID:${pair_id}\\tLB:Y\\tSM:${pair_id}\\tPL:ILLUMINA" $fasta $reads1 $reads2 \
    | elprep filter /dev/stdin ${pair_id}.bam \
    --mark-duplicates \
    --remove-duplicates \
    --filter-mapping-quality 0 \
    --clean-sam \
    --nr-of-threads ${task.cpus} \
    --sorting-order coordinate \
    --filter-unmapped-reads-strict
    samtools flagstats ${pair_id}.bam > ${pair_id}.bam_flagstats
    """
}

process variantcall {

    tag "$pair_id"
    publishDir params.outdir4, mode: 'copy'
    input:
    path fasta from params.fasta
    tuple val(pair_id), path(bamfile) from bam_ch
    output:
    set pair_id, "${pair_id}.vcf*"

    """
    samtools index ${pair_id}.bam
    freebayes -F 0.2 -u -p 2 -f $fasta $bamfile > ${pair_id}.vcf
    bgzip ${pair_id}.vcf
    tabix -p vcf ${pair_id}.vcf.gz
    """
}
