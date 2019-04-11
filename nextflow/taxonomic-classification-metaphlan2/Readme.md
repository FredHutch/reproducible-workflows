# Taxonomic classification of WGS reads with MetaPhlAn2

### Background 

MetaPhlAn2 is a great program for doing taxonomic classification of WGS microbiome datasets.


Truong DT, Franzosa EA, Tickle TL, Scholz M, Weingart G, Pasolli E, et al. MetaPhlAn2 for enhanced metagenomic taxonomic profiling. Nat Methods. 2015;12:902. https://doi.org/10.1038/nmeth.3589.

### Workflow

  1. Run MetaPhlan2 on every file in the manifest
  2. Write output to `output_folder`

### Input files

  * A reference database (by default included in the image)
  * A text file with one query FASTQ per line
