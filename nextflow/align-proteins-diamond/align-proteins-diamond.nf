#!/usr/bin/env nextflow

ref_fasta_f = file(params.ref_fasta)
query_ch = Channel.fromPath(file(params.sample_sheet).readLines())

process make_diamond_db {

  publishDir "$params.outdir/"

  container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"

  input:
  file fasta from ref_fasta_f

  output:
  file "${fasta.simpleName}.dmnd" into align_diamond

  """
  diamond makedb --in $fasta --db ${fasta.simpleName}.dmnd
  """
}

process align_diamond {

  publishDir "$params.outdir/"

  container "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"

  input:
  file ref from ref_fasta_f
  file query from query_ch

  output:
  file "${query.simpleName}.${ref.simpleName}.aln.gz"

  """
  set -e;

  gzip -t ${query} && \
    echo "Decompressing ${query}" && \
    gunzip -c ${query} > ${query}.fasta || \
    ln -s ${query} ${query}.fasta

  diamond \
      blastp \
      --db ${ref} \
      --query ${query}.fasta \
      --out ${query.simpleName}.${ref.simpleName}.aln \
      --outfmt 6 \
      --id 50 \
      --top 0 \
      --query-cover 50 \
      --subject-cover 50 \
      --threads 1;

  gzip ${query.simpleName}.${ref.simpleName}.aln
  """  
}
