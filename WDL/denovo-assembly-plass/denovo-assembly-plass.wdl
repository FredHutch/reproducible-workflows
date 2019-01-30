workflow denovoAssemblyPlass {

  File sample_sheet
  Array[File] sample_list = read_lines(sample_sheet)

  scatter (sample in sample_list) {
    call metaspades { 
      input: 
        input_fastq=sample
    }
  }

  output {
    Array[File] contigs=metaspades.contigs
  }

}

task metaspades {

  File input_fastq
  String memory="64G"
  String cpu="16"
  String output_name=basename(input_fastq, ".fastq.gz")
  
  runtime {
    docker: "quay.io/biocontainers/plass@sha256:c771c791ad89d9f2c09720d7e127d5b0e6ee2a35ca7688a1b79c461c116ddd05"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 

    plass assemble ${input_fastq} ${output_name}.faa tmp
    
    gzip ${output_name}.faa

  }
  output {
    File contigs = "${output_name}.faa.gz"
  }
}
