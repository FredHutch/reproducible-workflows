# Align Proteins with DIAMOND

### Background 

DIAMOND is a great drop-in replacement for BLASTp that achieves much better runtimes via double-indexing. 
[GH](https://github.com/bbuchfink/diamond)

Buchfink B, Xie C, Huson DH, "Fast and sensitive protein alignment using DIAMOND", Nature Methods 12, 59-60 (2015). 
[doi:10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)

### Workflow

  1. Use the reference FASTA to make a reference database
  2. Align each query FASTA against the reference database
  3. Combine all of the alignments into a single tarball

### Input files

  * One reference FASTA
  * A text file with one query FASTA per line

### Parameters
