# Megahit

### Background

_De novo_ assembly is the process of reconstructing complete and partial
genomes from whole-genome shotgun (WGS) sequence data. 
Megahit is a tool for performing _de novo_ assembly from WGS data.

### Workflow

#### Inputs

A sample sheet, listing each sample on a single line, with the location of
the FASTQ file for that sample. Each FASTQ file must be paired-end interleaved.

#### Outputs

An assembly for each individual sample, consisting of a FASTA file with the
contigs / scaffold for that assembly.

### Reference

Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015) MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, doi: 10.1093/bioinformatics/btv033 [PMID: 25609793].