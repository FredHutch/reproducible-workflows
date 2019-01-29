# MetaSPAdes

### Background

_De novo_ assembly is the process of reconstructing complete and partial
genomes from whole-genome shotgun (WGS) sequence data. 
MetaSPAdes is a tool for performing _de novo_ assembly from WGS data.

### Workflow

#### Inputs

A sample sheet, listing each sample on a single line, with the location of
the FASTQ file for that sample. Each FASTQ file must be paired-end interleaved.

#### Outputs

An assembly for each individual sample, consisting of a FASTA file with the
contigs / scaffold for that assembly.

### Reference
Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. Genome Res. 2017;27(5):824-834.
