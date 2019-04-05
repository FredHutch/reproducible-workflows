# Cluster Proteins with MMseqs2

### Background 

MMseqs2 is a great program for rapidly clustering proteins.


Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017).

Steinegger M and Soeding J. Clustering huge protein sequence sets in linear time. Nature Communications, doi: 10.1038/s41467-018-04964-5 (2018).

### Workflow

  1. Rename the FASTA headers in each file to include the sample name
  2. Concatenate the renamed FASTA files
  3. Cluster proteins with MMseqs2 and return the centroids and TSV with cluster information

### Input files

  * A set of reference FASTAs
  * A text file with one query FASTA per line

### Parameters

  * Clustering identity threshold
  * Clustering overlap threshold
