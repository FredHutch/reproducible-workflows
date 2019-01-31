workflow alignProteinsFAMLI {

  File sample_sheet
  File ref_fastp
  Array[File] sample_list = read_lines(sample_sheet)

  call MakeDiamondDatabase {
    input:
      fasta=ref_fastp
  }

  scatter (sample in sample_list) {
    call DiamondBlastx { 
      input: 
        input_fastq=sample,
        refdb=MakeDiamondDatabase.db
    }
    call FAMLI { 
      input: 
        input_aln=DiamondBlastx.aln
    }
  }

  output {
    Array[File] results=FAMLI.results
  }

}

task MakeDiamondDatabase {

  File fasta
  String fasta_base = basename(fasta, ".faa.gz")

  runtime {
    docker: "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    memory: "2G"
    cpu: "1"
  }

  command {
    diamond makedb --in ${fasta} --db ${fasta_base}.db.dmnd
  }
  output {
    File db = "${fasta_base}.db.dmnd"
  }
}

task DiamondBlastx {

  File refdb
  File input_fastq
  String db_prefix = sub(refdb, ".dmnd", "")
  String query_name = basename(input_fastq, ".fastq.gz")
  String min_id = "80"
  String top_pct = "10"
  String query_cover = "95"
  String subject_cover = "50"
  String cpu = "1"
  String min_score="20"
  String blocks="15"
  String query_gencode="11"
  
  runtime {
    docker: "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    memory: "16G"
    cpu: cpu
  }

  command {
    set -e;

    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${query_name}.aln \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover} \
      --id ${min_id} \
      --top ${top_pct} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
      --unal 0
  }

  output {
    File aln = "${query_name}.aln"
  }
}

task FAMLI {

  File input_aln
  String memory="16G"
  String cpu="16"
  String output_name=basename(input_aln, ".aln")
  String batchsize="50000000"
  
  runtime {
    docker: "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 
    
    famli \
      filter \
      --input ${input_aln} \
      --output ${output_name}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}

    gzip ${output_name}.json

  }
  output {
    File results = "${output_name}.json.gz"
  }
}
