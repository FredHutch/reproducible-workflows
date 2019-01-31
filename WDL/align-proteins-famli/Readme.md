# FAMLI

### Background

FAMLI aligns a set of reads (FASTQ) against a reference database of protein-coding genes,
filters the resulting alignment, and reports the relative abundance of each gene.

### Workflow

#### Inputs

A sample sheet, listing each sample on a single line, with the location of
the FASTQ file for that sample. Each FASTQ file must be paired-end interleaved.
Reference database must be indexed for DIAMOND.

#### Outputs

A JSON with the abundance of each gene
