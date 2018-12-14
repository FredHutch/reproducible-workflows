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

  output {
    File integrated_assembly_hdf5=IntegrateAssemblies.integrated_assembly_hdf5
    File integrated_assembly_dmnd=IntegrateAssemblies.integrated_assembly_dmnd
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
    spades.py \
        -s ${input_file} \
        -o ${sample_name} \
        -t 16 \
        --phred-offset 33 \
        -m 64; 
    mv output/scaffolds.fasta ${sample_name}.scaffolds.fasta
    gzip ${sample_name}.scaffolds.fasta
  }
  output {
    File assembly = "${sample_name}.scaffolds.fasta.gz"
  }
}

task Annotation {

  File assembly

  runtime {
    docker: "quay.io/fhcrc-microbiome/metaspades:v3.11.1--8"
    memory: "16G"
    cpu: "4"
  }

  command {
    sample_name="$(echo ${assembly} | sed 's/.*\///' | sed 's/.scaffolds.fasta.gz//')"
    gunzip -c "${assembly}" > scaffolds.fasta; 
    prokka \
        --outdir ./ \
        --prefix "$sample_name" \
        --cpus 4 \
        --metagenome \
        scaffolds.fasta; 
    gzip "$sample_name.fastp";
    gzip "$sample_name.gff"  }
  output {
    File fastp = "$sample_name.fastp"
    File gff = "$sample_name.gff"
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
    mkdir assembly_files;
    for f in ${sep=" " fastp_list}; do cp "$f" assembly_files/; done; 
    for f in ${sep=" " gff_list}; do cp "$f" assembly_files/; done; 
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
