# Plass

### Background

_De novo_ assembly is the process of reconstructing complete and partial
genomes from whole-genome shotgun (WGS) sequence data. 
Plass is a tool for performing _de novo_ assembly of protein-coding genes from WGS data.

### Workflow

#### Inputs

A sample sheet, listing each sample on a single line, with the location of
the FASTQ file for that sample. Each FASTQ file must be paired-end interleaved.

#### Outputs

An assembly for each individual sample, consisting of a FASTP file with the
protein-coding genes for that assembly.

### Reference
Protein-level assembly increases protein sequence recovery from metagenomic samples manyfold
Martin Steinegger, Milot Mirdita, Johannes Soding
bioRxiv 386110; doi: https://doi.org/10.1101/386110
This article is a preprint and has not been peer-reviewed
