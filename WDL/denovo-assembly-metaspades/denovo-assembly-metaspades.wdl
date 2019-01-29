workflow denovoAssemblyMetaspades {

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
    docker: "quay.io/biocontainers/spades@sha256:9f097c5d6d7944b68828e10d94504ac49a93bf337a9afed17232594b126b807e"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 

    spades_memory=$(echo ${memory} | tr -d 'G')
    
    spades.py \
        --12 ${input_fastq} \
        -o ${output_name} \
        --threads ${cpu} \
        --memory $spades_memory \
        --meta \
        --phred-offset 33 \
        --only-assembler;
    
    mv ${output_name}/scaffolds.fasta ${output_name}.scaffolds.fasta
    gzip ${output_name}.scaffolds.fasta

  }
  output {
    File contigs = "${output_name}.scaffolds.fasta.gz"
  }
}
