#!/bin/bash

set -e

nextflow \
    run \
    align-proteins-diamond.nf \
    --sample_sheet local_inputs/sample_sheet.txt \
    --ref_fasta local_inputs/ncbi_ref_proteins.faa \
    --outdir local_outputs \
    -with-docker "ubuntu:16.04" \
    -resume
