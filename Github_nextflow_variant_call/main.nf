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
println "genome: $params.genome"

// process

/*
 * Step 1 cleaning
 */
process preprocessing {
    cpus 2
    tag "$pair_id"
    publishDir params.outdir1, mode: 'copy'
    input:
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    set pair_id, "${pair_id}_QT_R1.fq.gz" into read1_ch
    set pair_id, "${pair_id}_QT_R2.fq.gz" into read2_ch
    set pair_id, "${pair_id}_fastp*"

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
process breseq {

    tag "$pair_id"
    publishDir params.outdir2, mode: 'copy'
    cpus 4
    input:
    path gff from params.gff
    tuple val(pair_id), path(reads1) from read1_ch
    tuple val(pair_id), path(reads2) from read2_ch

    output:
    set pair_id, "${pair_id}"


    """
    breseq -j ${task.cpus} -o ${pair_id} -r $gff $reads1 $reads2
    rm -rf ${pair_id}/0* ${pair_id}/data
    mv ${pair_id}/output/* ${pair_id}/
    mv ${pair_id}/output.gd ${pair_id}/${pair_id}.gd
    mv ${pair_id}/output.vcf ${pair_id}/${pair_id}.vcf
    mv ${pair_id}/summary.html ${pair_id}/${pair_id}.html
    mv ${pair_id}/summary.json ${pair_id}/${pair_id}.json
    """
}
