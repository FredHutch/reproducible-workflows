workflow denovoAssemblyMegahit {

  File sample_sheet
  Array[File] sample_list = read_lines(sample_sheet)

  scatter (sample in sample_list) {
    call megahit { 
      input: 
        input_fastq=sample
    }
  }

  output {
    Array[File] contigs=megahit.contigs
  }

}

task megahit {

  File input_fastq
  String memory="64G"
  String cpu="16"
  String output_name=basename(input_fastq, ".fastq.gz")
  
  runtime {
    docker: "quay.io/biocontainers/megahit@sha256:8c9f17dd0fb144254e4d6a2a11d46b522239d752d2bd15ae3053bb1a31cc6d01"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 
    
    megahit --12 ${input_fastq} \
            -o ${output_name} \
            -t ${cpu};

    mv ${output_name}/final.contigs.fa ${output_name}.contigs.fasta
    gzip ${output_name}.contigs.fasta

  }
  output {
    File contigs = "${output_name}.contigs.fasta.gz"
  }
}
