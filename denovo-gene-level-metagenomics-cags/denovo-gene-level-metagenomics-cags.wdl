workflow denovoGeneLevelMetagenomicsCAGs {

  File manifest
  Array[Array[String]] inputSamples = read_tsv(manifest)

  scatter (sample in inputSamples) {
    call deNovoAssembly { 
      input: input_file=sample[1], 
             sample_name=sample[0]
    }
  }
  scatter (assembly in deNovoAssembly.assembly) {
    call Annotation { 
      input: assembly=assembly 
    }
  }

  call IntegrateAssemblies { 
    input: fastp_list=Annotation.fastp, 
        gff_list=Annotation.gff
  }

  scatter (sample in inputSamples) {
    call AlignDiamond { 
      input: input_file=sample[1], 
             sample_name=sample[0],
             ref=IntegrateAssemblies.integrated_assembly_dmnd
    }
  }

  scatter (alignment in AlignDiamond.alignment) {
    call FilterFAMLI { 
      input: alignment=alignment
    }
  }

  output {
    File integrated_assembly_hdf5=IntegrateAssemblies.integrated_assembly_hdf5
    File integrated_assembly_dmnd=IntegrateAssemblies.integrated_assembly_dmnd
    Array[File] famli_abundance=FilterFAMLI.abundance
  }

}

task deNovoAssembly {

  File input_file
  String sample_name

  runtime {
    docker: "quay.io/fhcrc-microbiome/metaspades:v3.11.1--10"
    memory: "64G"
    cpu: "16"
  }

  command {
    set -e; 
    spades.py \
        -s ${input_file} \
        -o ${sample_name} \
        -t 16 \
        --phred-offset 33 \
        -m 64; 
    mv ${sample_name}/scaffolds.fasta ${sample_name}.scaffolds.fasta
    gzip ${sample_name}.scaffolds.fasta
  }
  output {
    File assembly = "${sample_name}.scaffolds.fasta.gz"
  }
}

task Annotation {

  File assembly
  String sample_name = sub(sub(assembly, ".*/", ""), ".scaffolds.fasta.gz$", "")

  runtime {
    docker: "quay.io/fhcrc-microbiome/metaspades:v3.11.1--8"
    memory: "16G"
    cpu: "4"
  }

  command {
    set -e; 
    gunzip -c "${assembly}" > scaffolds.fasta; 
    prokka \
        --outdir output/ \
        --prefix "${sample_name}" \
        --cpus 4 \
        --metagenome \
        scaffolds.fasta; 
    mv "output/${sample_name}.faa" "output/${sample_name}.fastp";
    gzip "output/${sample_name}.fastp";
    gzip "output/${sample_name}.gff"  
  }
  output {
    File fastp = "output/${sample_name}.fastp.gz"
    File gff = "output/${sample_name}.gff.gz"
  }
}

task IntegrateAssemblies {

  Array[File] fastp_list
  Array[File] gff_list

  runtime {
    docker: "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    memory: "16G"
    cpu: "4"
  }

  command {
    set -e; 
    mkdir assembly_files;
    echo ${sep=" " fastp_list}; 
    echo ${sep=" " gff_list}; 
    for f in ${sep=" " fastp_list}; do ln -s "$f" assembly_files/; done; 
    for f in ${sep=" " gff_list}; do ln -s "$f" assembly_files/; done; 
    integrate_assemblies.py \
        --gff-folder assembly_files/ \
        --prot-folder assembly_files/ \
        --output-name integrated_assembly \
        --output-folder ./ \
        --gff-suffix gff \
        --prot-suffix fastp \
        --temp-folder ./
  }
  output {
    File integrated_assembly_hdf5 = "integrated_assembly.hdf5"
    File integrated_assembly_dmnd = "integrated_assembly.dmnd"
  }
}

task AlignDiamond {

  File input_file
  String sample_name
  File ref

  runtime {
    docker: "quay.io/fhcrc-microbiome/famli:v1.1"
    memory: "16G"
    cpu: "4"
  }

  command {
    set -e; 
    diamond \
        blastx \
        --query "${input_file}" \
        --out "${sample_name}.aln" \
        --threads 4 \
        --db  ${ref} \
        --outfmt  6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --min-score  20 \
        --query-cover  95 \
        --id 80 \
        --top 10 \
        --block-size  4 \
        --query-gencode 11 \
        --unal 0; 
    gzip "${sample_name}.aln"
  }
  output {
    File alignment = "${sample_name}.aln.gz"
  }
}

task FilterFAMLI {

  File alignment
  String output_fp = sub(sub(alignment, ".*/", ""), ".aln.gz", ".json")

  runtime {
    docker: "quay.io/fhcrc-microbiome/famli:v1.1"
    memory: "16G"
    cpu: "4"
  }

  command {
    set -e; 
    famli \
        filter \
        --input "${alignment}" \
        --output "${output_fp}" \
        --threads 4 \
        --batchsize 50000000

    gzip "${output_fp}"
  }
  output {
    File abundance = "${output_fp}.gz"
  }
}

