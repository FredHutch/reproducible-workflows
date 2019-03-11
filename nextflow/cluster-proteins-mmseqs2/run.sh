#!/bin/bash

set -e

nextflow \
    run \
    cluster-proteins-mmseqs2.nf \
    --sample_sheet local_inputs/sample_sheet.txt \
    --outdir local_outputs \
    -with-docker "ubuntu:16.04" \
    -resume
