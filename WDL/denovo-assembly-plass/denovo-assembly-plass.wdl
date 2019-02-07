workflow denovoAssemblyPlass {

  File sample_sheet
  String translation_table = "11"
  Array[File] sample_list = read_lines(sample_sheet)

  scatter (sample in sample_list) {
    call plass { 
      input: 
        input_fastq=sample,
        translation_table=translation_table
    }
  }

  output {
    Array[File] contigs=plass.contigs
  }

}

task plass {

  File input_fastq
  String translation_table = "11"
  String min_orf_length = "20"
  String memory="64G"
  String cpu="16"
  String output_name=basename(input_fastq, ".fastq.gz")
  
  runtime {
    # docker: "quay.io/biocontainers/plass@sha256:c771c791ad89d9f2c09720d7e127d5b0e6ee2a35ca7688a1b79c461c116ddd05"
    docker: "soedinglab/plass"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 

    plass assemble --use-all-table-starts --min-length "${min_orf_length}" --translation-table ${translation_table} "${input_fastq}" "${output_name}.faa" tmp

    rm -r tmp
    
    gzip ${output_name}.faa

  }
  output {
    File contigs = "${output_name}.faa.gz"
  }
}
