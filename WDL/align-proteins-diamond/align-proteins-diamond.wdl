workflow alignProteinsDiamond {

  File query_list
  File reference
  Array[File] queryList = read_lines(query_list)

  call MakeDiamondDatabase {
      input: fasta=reference
  }

  scatter (query in queryList) {
    call RunDiamond { 
      input: 
        db=MakeDiamondDatabase.db, 
        query=query
    }
  }

  call MakeTarball { 
    input: files=RunDiamond.aln
  }

  output {
    File combined_alignments=MakeTarball.combined_alignments
  }

}

task MakeDiamondDatabase {

  File fasta
  String fasta_base = sub(fasta, ".*/", "")

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

task RunDiamond {

  File db
  File query
  String db_prefix = sub(db, ".dmnd", "")
  String db_name = sub(sub(db, ".dmnd", ""), ".*/", "")
  String query_name = sub(query, ".*/", "")

  runtime {
    docker: "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    memory: "2G"
    cpu: "1"
  }

  command {
    set -e;

    diamond \
        blastp \
        --db ${db} \
        --query ${query} \
        --out ${query_name}.${db_name}.aln \
        --outfmt 6 \
        --id 50 \
        --max-target-seqs 500 \
        --query-cover 50 \
        --subject-cover 50 \
        --threads 1;

    gzip ${query_name}.${db_name}.aln
  }
  output {
    File aln = "${query_name}.${db_name}.aln.gz"
  }
}

task MakeTarball {

  Array[File] files

  runtime {
    docker: "quay.io/fhcrc-microbiome/docker-diamond@sha256:0f06003c4190e5a1bf73d806146c1b0a3b0d3276d718a50e920670cf1bb395ed"
    memory: "2G"
    cpu: "1"
  }

  command {
    set -e; 
    tar cvf combined_alignments.tar ${sep=" " files};
  }
  output {
    File combined_alignments = "combined_alignments.tar"
  }
}
