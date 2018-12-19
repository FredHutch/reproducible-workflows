workflow denovoGeneLevelMetagenomicsCAGs {

  File manifest
  File tax_db
  File taxonmap
  File taxonnodes
  Array[Array[String]] inputSamples = read_tsv(manifest)
  String cags_min_samples = "2"
  String cags_normalization = "sum"
  String cags_max_dist = "0.2"

  scatter (sample in inputSamples) {
    call deNovoAssembly { 
      input: 
        input_file=sample[1], 
        sample_name=sample[0]
    }

    call Annotation { 
      input: assembly=deNovoAssembly.assembly
    }
  }

  call IntegrateAssemblies { 
    input: fastp_list=Annotation.fastp, 
        gff_list=Annotation.gff
  }

  call AnnotateProteinsTax {
    input: 
      ref=tax_db,
      taxonmap=taxonmap,
      taxonnodes=taxonnodes,
      query=IntegrateAssemblies.integrated_assembly_fastp
  }

  scatter (sample in inputSamples) {
    call AlignDiamond { 
      input: 
        input_file=sample[1], 
        sample_name=sample[0],
        ref=IntegrateAssemblies.integrated_assembly_dmnd
    }

    call FilterFAMLI { 
      input: alignment=AlignDiamond.alignment
    }
  }

  call MakeCAGs {
    input: 
      filtered_gene_abundances=FilterFAMLI.abundance,
      cags_min_samples=cags_min_samples,
      cags_normalization=cags_normalization,
      cags_max_dist=cags_max_dist
  }

  output {
    File integrated_assembly_hdf5=IntegrateAssemblies.integrated_assembly_hdf5
    File integrated_assembly_dmnd=IntegrateAssemblies.integrated_assembly_dmnd
    Array[File] famli_abundance=FilterFAMLI.abundance
    File integrated_assembly_tax=AnnotateProteinsTax.tax_annot
    File cags_json=MakeCAGs.cags_json
  }

}

task deNovoAssembly {

  File input_file
  String sample_name

  runtime {
    docker: "quay.io/fhcrc-microbiome/metaspades@sha256:5cf9216f3536c1b1f1f214972ba89d358033cfea5d28ce90b8a3fb7764eccda5"
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
    docker: "quay.io/fhcrc-microbiome/metaspades@sha256:cb5885c0f37bb2a264c81753fb060212db98c39d9962464616a395aab5fdde15"
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
    docker: "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies@sha256:dc8b0ffc62126f27219708e4b311f23878faca9e4a462e3c8ef07632ed11e531"
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
    File integrated_assembly_fastp = "integrated_assembly.fastp.gz"
    File integrated_assembly_hdf5 = "integrated_assembly.hdf5"
    File integrated_assembly_dmnd = "integrated_assembly.dmnd"
  }
}

task AnnotateProteinsTax {

  File query
  File ref
  File taxonmap
  File taxonnodes
  String query_base = sub(sub(query, ".*/", ""), ".fastp.gz", "")
  String ref_base = sub(sub(ref, ".*/", ""), ".dmnd", "")
  String top_pct = "1"
  String blocks = "5"
  String threads = "4"

  runtime {
    docker: "quay.io/fhcrc-microbiome/famli@sha256:cd57ca9edd302e3b0803b8374e6e1d70eaa32edff9210f8b6dff28d373835b22"
    memory: "16G"
    cpu: "4"
  }

  command {
    set -e; 
    diamond \
      blastp \
      --db ${ref} \
      --query ${query} \
      --out ${query_base}.${ref_base}.tax \
      --taxonmap ${taxonmap} \
      --taxonnodes ${taxonnodes} \
      --outfmt 102 \
      --top ${top_pct} \
      -b ${blocks} \
      --threads ${threads}
    gzip "${query_base}.${ref_base}.tax"
  }
  output {
    File tax_annot = "${query_base}.${ref_base}.tax.gz"
  }
}

task AlignDiamond {

  File input_file
  String sample_name
  File ref

  runtime {
    docker: "quay.io/fhcrc-microbiome/famli@sha256:cd57ca9edd302e3b0803b8374e6e1d70eaa32edff9210f8b6dff28d373835b22"
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
    docker: "quay.io/fhcrc-microbiome/famli@sha256:cd57ca9edd302e3b0803b8374e6e1d70eaa32edff9210f8b6dff28d373835b22"
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

task MakeCAGs {
  Array[File] filtered_gene_abundances
  String cags_min_samples = "2"
  String cags_normalization = "sum"
  String cags_max_dist = "0.2"

  runtime {
    docker: "quay.io/fhcrc-microbiome/find-cags@sha256:bfde01058cb80e970818853a8d4f64eb3184138e84671fc12819f348a38deedf"
    memory: "120G"
    cpu: "16"
  }

  command {
    set -e;

    python -c "
import json;
json.dump(dict([
  (f, f)
  for f in '${sep=' ' filtered_gene_abundances}'.split(' ')
]), open('sample_sheet.json', 'wt'))"

    cat sample_sheet.json; 

    for f in ${sep=" " filtered_gene_abundances}; do
      echo $f;
      gunzip -c $f
    done

    find-cags.py \
      --sample-sheet sample_sheet.json \
      --output-prefix output \
      --output-folder ./ \
      --normalization ${cags_normalization} \
      --max-dist ${cags_max_dist} \
      --min-samples ${cags_min_samples}

  }
  output {
    File cags_json = "output.cags.json.gz"
  }

}

