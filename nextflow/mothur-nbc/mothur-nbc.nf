#!/usr/bin/env nextflow

Channel.fromPath("${params.input_directory}*.gz")
       .set{ fasta_ch }

process mothurNBC {

  container "quay.io/fhcrc-microbiome/mothur@sha256:02c9eba4525267343a9b82eb1b6490cdf112d5d016867568e3f38a23425bc5ea"
  cpus 4
  memory "16 GB"
  publishDir "${params.output_directory}"
  errorStrategy 'retry'

  input:
  file fasta from fasta_ch

  output:
  file "${fasta}.tax.counts"

  """
set -e

gunzip -c ${fasta} > input.fasta
gunzip -k -f /usr/local/dbs/silva.bacteria.fasta.gz

echo "classify.seqs(fasta=input.fasta, template=/usr/local/dbs/silva.bacteria.fasta, taxonomy=/usr/local/dbs/silva.bacteria.gg.tax, method=wang, ksize=8, iters=100, processors=4)" > batchfile

mothur batchfile

cat *.wang.taxonomy | cut -f 2 | sort | uniq -c | sort -nr > ${fasta}.tax.counts

  """
}
