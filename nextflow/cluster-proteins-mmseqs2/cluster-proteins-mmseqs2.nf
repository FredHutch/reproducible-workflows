#!/usr/bin/env nextflow

fasta_ch = Channel.fromPath(file(params.sample_sheet).readLines())
params.identity = "0.90"
params.overlap = "0.50"

process rename_fasta_headers {

  container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
  input:
  file fasta from fasta_ch
  output:
  file "${fasta}.renamed.fasta" into renamed_fasta_ch

  """
  set -e;

python3 << END
from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip

def gzip_safe_open(fp, mode="rt"):
    if fp.endswith(".gz"):
        return gzip.open(fp, mode)
    else:
        return open(fp, mode)

with open("${fasta}.renamed.fasta", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip_safe_open("${fasta}")):
        header = header.split(" ")[0]
        fo.write(">" + "${fasta}" + ":" + header + "\\n" + seq + "\\n")
END
  """
}

process concatenate_fastas {
  container "ubuntu:16.04"
  input:
  file fasta_list from renamed_fasta_ch.collect()
  output:
  file "all.fasta.gz" into concat_fasta_ch

  """
  for fasta in ${fasta_list}
    do
        cat \$fasta >> ALL_FASTAS
    done

  mv ALL_FASTAS all.fasta
  gzip all.fasta
  """
}

process mmseqs2 {
    container "quay.io/biocontainers/mmseqs2@sha256:f935cdf9a310118ba72ceadd089d2262fc9da76954ebc63bafa3911240f91e06"
    input:
    file fasta_in from concat_fasta_ch
    val identity from params.identity
    val overlap from params.overlap
    output:
    file "${fasta_in}.rep.fasta"
    file "${fasta_in}.clusters.tsv"
    publishDir params.outdir

    """
    set -e;

    mmseqs createdb ${fasta_in} DB;
    mmseqs cluster DB clu tmp --min-seq-id ${identity} --max-seqs 1000000 -c ${overlap} --threads 1;
    mmseqs createtsv DB DB clu ${fasta_in}.clusters.tsv;
    mmseqs result2repseq DB clu ${fasta_in}.rep;
    mmseqs result2flat DB DB ${fasta_in}.rep ${fasta_in}.rep.fasta --use-fasta-header

    # Make sure the output file has data
    echo "Number of input sequences: \$(gunzip -c ${fasta_in} | grep -c '>')"
    echo "Number of clustered sequences: \$(cat "${fasta_in}.clusters.tsv" | wc -l)"
    echo "Number of clusters: \$(cat "${fasta_in}.clusters.tsv" | cut -f 1 | sort -u | wc -l)"
    echo "Number of representative sequences: \$(cat "${fasta_in}.rep.fasta" | grep -c '>')"
    (( \$(cat "${fasta_in}.rep.fasta" | grep -c '>') == \$(cat "${fasta_in}.clusters.tsv" | cut -f 1 | sort -u | wc -l) ))
    (( \$(cat "${fasta_in}.clusters.tsv" | wc -l) == \$(gunzip -c ${fasta_in} | grep -c '>') ))

    """

}