workflow genomeAnnotationProkka {

  File sample_sheet
  Array[File] sample_list = read_lines(sample_sheet)

  scatter (sample in sample_list) {
    call prokka { 
      input: 
        input_fasta=sample
    }
  }

  output {
    Array[File] gff=prokka.gff
    Array[File] faa=prokka.faa
    Array[File] gbk=prokka.gbk
  }

}

task prokka {

  File input_fasta
  String memory="16G"
  String cpu="16"
  String output_name=basename(input_fasta, ".fasta.gz")
  
  runtime {
    docker: "quay.io/biocontainers/prokka@sha256:6005120724868b80fff0acef388de8c9bfad4917b8817f383703eeacd979aa5a"
    memory: memory
    cpu: cpu
  }

  command {
    set -e; 
    
    gunzip -c ${input_fasta} > input.fasta; 
    prokka --outdir ${output_name} --prefix ${output_name} --metagenome input.fasta
    gzip ${output_name}/${output_name}.*

  }
  output {
    File gff = "${output_name}/${output_name}.gff.gz"
    File faa = "${output_name}/${output_name}.faa.gz"
    File gbk = "${output_name}/${output_name}.gbk.gz"
  }
}
