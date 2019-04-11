#!/bin/bash

set -e

nextflow \
    run \
    metaphlan2.nf \
    --manifest local_inputs/sample_sheet.txt \
    --output_folder local_outputs \
    -with-docker "ubuntu:16.04" \
    -resume
