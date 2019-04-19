#!/bin/bash

set -e

nextflow \
    run \
    interleave-fastq-pairs.nf \
    --input_folder s3://fh-ctr-public-reference-data/workflow_testing_data/NF/interleave-fastq-pairs/input/ \
    --output_folder s3://fh-ctr-public-reference-data/workflow_testing_data/NF/interleave-fastq-pairs/output/ \
    -resume
