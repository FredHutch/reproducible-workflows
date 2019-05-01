#!/usr/bin/env nextflow

fastq_ch = Channel.from(file(params.manifest).readLines())
                  .map {it -> file(it)}
params.input_type = "fastq"

process metaphlan2 {
    container "quay.io/fhcrc-microbiome/metaphlan@sha256:51b416458088e83d0bd8d840a5a74fb75066b2435d189c5e9036277d2409d7ea"
    cpus 16
    memory "32 GB"
    publishDir "${params.output_folder}"

    input:
    file input_fastq from fastq_ch
    val input_type from params.input_type

    output:
    file "${input_fastq}.metaphlan.tsv"

    """
    metaphlan2.py --input_type ${input_type} --tmp_dir ./ -o ${input_fastq}.metaphlan.tsv ${input_fastq}
    """
}