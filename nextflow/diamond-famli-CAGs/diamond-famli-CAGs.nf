#!/usr/bin/env nextflow

Channel
    .fromFilePairs("${params.input_folder}*_R{1,2}_*.fastq.gz")
    .ifEmpty { error "Cannot find any reads matching: ${params.input_folder}*_R{1,2}_*.fastq.gz" }
    .set { fastq_pair_ch }

process interleave {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "120 GB"
    errorStrategy "retry"
    publishDir params.output_folder
    
    input:
    set prefix, file(reads) from fastq_pair_ch
    
    output:
    file "${prefix}.fastq.gz" into diamond_fastq_ch, metaphlan_fastq_ch

"""
#!/usr/bin/env python
import gzip
import os

read1, read2 = "${reads}".split(" ", 1)

assert os.path.exists(read1), [read1, os.listdir(".")]
assert os.path.exists(read2), [read2, os.listdir(".")]

def split_whitespace(line, suffix=None):
    line = line.rstrip("\\n").split("\\t")[0].split(" ")[0]
    if suffix is not None:
        return line + suffix
    return line

n_reads = 0
with gzip.open(read1, "rt") as f1, gzip.open(read2, "rt") as f2, gzip.open("${prefix}.fastq.gz", "wt") as fo:
    while True:
        line = f1.readline()
        if line.strip() == "":
            break

        fo.write(split_whitespace(line, suffix="/1\\n"))
        
        for i in range(3):
            fo.write(f1.readline())
        
        fo.write(split_whitespace(f2.readline(), suffix="/2\\n"))
        for i in range(3):
            fo.write(f2.readline())

        n_reads += 1
        if n_reads % 1000000 == 0:
            print("Processed " + str(n_reads) + " paired reads")

"""    
}


process metaphlan2 {
    container "quay.io/fhcrc-microbiome/metaphlan@sha256:51b416458088e83d0bd8d840a5a74fb75066b2435d189c5e9036277d2409d7ea"
    cpus 16
    memory "32 GB"
    publishDir "${params.output_folder}"

    input:
    file input_fastq from metaphlan_fastq_ch
    
    output:
    file "${input_fastq}.metaphlan.tsv"

    """
    metaphlan2.py --input_type fastq --tmp_dir ./ -o ${input_fastq}.metaphlan.tsv ${input_fastq}
    """
}


process diamond {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    cpus 32
    memory "240 GB"
    errorStrategy "retry"
    
    input:
    file refdb from file(params.diamond_db)
    file input_fastq from diamond_fastq_ch
    val min_id from 90
    val query_cover from 50
    val cpu from 32
    val top from 1
    val min_score from 20
    val blocks from 15
    val query_gencode from 11

    output:
    file "${input_fastq}.aln.gz" into aln_ch

    """
    set -e;

    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${input_fastq}.aln.gz \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover} \
      --id ${min_id} \
      --top ${top} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
      --compress 1 \
      --unal 0

    # Removing reference database and input FASTQ
    rm ${refdb} ${input_fastq}
    
    echo "Done"
    """

}


process famli {
    container "quay.io/fhcrc-microbiome/famli@sha256:241a7db60cb735abd59f4829e8ddda0451622b6eb2321f176fd9d76297d8c9e7"
    cpus 16
    memory "120 GB"
    errorStrategy "retry"
    publishDir "${params.output_folder}"
    
    input:
    file input_aln from aln_ch
    val cpu from 16
    val batchsize from 50000000

    output:
    file "${input_aln.simpleName}.json.gz" into famli_json_for_filter, famli_json_for_summary

    """
    set -e; 
    
    famli \
      filter \
      --input ${input_aln} \
      --output ${input_aln.simpleName}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}

    gzip ${input_aln.simpleName}.json
    """

}


process extract_reference_sequences {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    cpus 1
    memory "4 GB"
    
    input:
    file refdb from file(params.diamond_db)

    output:
    file "refdb.fasta.gz" into refdb_fasta

    """
    diamond getseq -d ${refdb} > refdb.fasta
    gzip refdb.fasta
    """

}


process filter_detected_genes {
    container "quay.io/biocontainers/biopython@sha256:1196016b05927094af161ccf2cd8371aafc2e3a8daa51c51ff023f5eb45a820f"
    cpus 16
    memory "32 GB"
    publishDir params.output_folder
    
    input:
    file refdb_fasta
    file famli_json_list from famli_json_for_filter.collect()
    
    output:
    file "refdb.filtered.fasta.gz" into refdb_filtered_fasta

    """
#!/usr/bin/env python

import os
import json
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Make a set of all of the genes that were detected
detected_genes = set([])
for f in os.listdir("."):
    if f.endswith(".json.gz"):
        print(f)
        for r in json.load(gzip.open(f, "rt")):
            detected_genes.add(r["id"].split(" ")[0])
print("Detected " + str(len(detected_genes)) + " genes")

n = 0
with gzip.open("refdb.filtered.fasta.gz", "wt") as fo:
    for header, seq in SimpleFastaParser(gzip.open("${refdb_fasta}", "rt")):
        header = header.split(" ")[0]
        if header in detected_genes:
            n += 1
            fo.write(">" + header + "\\n" + seq + "\\n")

print("Wrote out " + str(n) + " detected genes")
assert n == len(detected_genes)

    """

}

process taxonomic_annotation {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    cpus 32
    memory "240 GB"
    publishDir params.output_folder

    input:
    file query from refdb_filtered_fasta
    file diamond_tax_db from file(params.diamond_tax_db)
    val min_id from 90
    val blocks from 20
    val cpu from 16
    
    output:
    file "refdb.tax.aln.gz" into refdb_tax

    """
    diamond \
      blastp \
      --db ${diamond_tax_db} \
      --query ${query} \
      --out refdb.tax.aln.gz \
      --outfmt 102 \
      --id ${min_id} \
      --top ${100 - min_id} \
      --block-size ${blocks} \
      --threads ${cpu} \
      --compress 1

    rm ${diamond_tax_db}
    """

}


process eggnog_annotation {
    container "quay.io/biocontainers/eggnog-mapper:1.0.3--py27_0"
    cpus 16
    memory "120 GB"
    publishDir params.output_folder
    
    input:
    file query from refdb_filtered_fasta
    file db from file(params.eggnog_db)
    file dmnd_db from file(params.eggnog_dmnd_db)

    output:
    file "${query.simpleName}.emapper.annotations.gz" into refdb_eggnog
    
    """
    mkdir data
    mkdir TEMP
    mkdir SCRATCH
    mv ${db} data/eggnog.db
    mv ${dmnd_db} data/eggnog_proteins.dmnd

    emapper.py \
      -i ${query} \
      --output ${query.simpleName} \
      -m "diamond" \
      --cpu 16 \
      --data_dir data/ \
      --scratch_dir SCRATCH/ \
      --temp_dir TEMP/ \

    gzip ${query.simpleName}.emapper.annotations
    
    """

}
