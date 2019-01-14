workflow Panel_BWA_GATK4_Samtools_Var_Annotate {
  File batchFile
  Array[Object] batchInfo = read_objects(batchFile)

  String ref_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices


scatter (job in batchInfo){
  String sampleName = job.sampleName
  File bamLocation = job.bamLocation
  File bedLocation = job.bedLocation

# Get the basename, i.e. strip the filepath and the extension
  String bam_basename = basename(bamLocation, ".bam")
  String base_file_name = sampleName + "." + ref_name

# convert unmapped bam to fastq
  call SamToFastq {
    input:
      input_bam = bamLocation,
      base_file_name = base_file_name,
  }
  # Map reads to reference
  call BwaMem {
    input:
      input_fastq = SamToFastq.output_fastq,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa
  }
  # # Merge original uBAM and BWA-aligned BAM
  # call MergeBamAlignment {
  #   input:
  #     unmapped_bam = bamLocation,
  #     bwa_commandline = bwa_commandline,
  #     bwa_version = BWAVersion,
  #     aligned_bam = SamToFastqAndBwaMem.output_bam,
  #     output_bam_basename = bam_basename + ".aligned.unsorted",
  #     ref_fasta = ref_fasta,
  #     ref_fasta_index = ref_fasta_index,
  #     ref_dict = ref_dict,
  #     ref_alt = ref_alt,
  #     ref_amb = ref_amb,
  #     ref_ann = ref_ann,
  #     ref_bwt = ref_bwt,
  #     ref_pac = ref_pac,
  #     ref_sa = ref_sa

  # }
  # mark duplicates
  # call MarkDuplicates {
  #   input:
  #     input_bams = MergeBamAlignment.output_bam,
  #     output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
  #     metrics_filename = base_file_name + ".duplicate_metrics",
  #     docker_image = gatk_docker,
  #     gatk_path = gatk_path,
  #     disk_size = large_disk,
  #     compression_level = compression_level,
  #     preemptible_tries = preemptible_tries
  # }



 # End scatter 
}
  # End workflow
}

# TASK DEFINITIONS

# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  File input_bam
  String base_file_name

  command {
    set -e

    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=${base_file_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true \

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14G"
    cpu: "1"
  }
  output {
    File output_fastq = "${base_file_name}.fastq"
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task BwaMem {
  File input_fastq
  String base_file_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  command {
    set -o pipefail
    set -e

    /usr/gitc/bwa mem \
      -p -v 3 -t 16 \
      -M \
      ${ref_fasta} ${input_fastq} | samtools view \
      -@16 -b -1 > ${base_file_name}.bam \
      2> ${base_file_name}.bwa.stderr.log

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14G"
    cpu: "16"
  }
  output {
    File output_bam = "${base_file_name}.unmapped.bam"
    File bwa_stderr_log = "${base_file_name}.bwa.stderr.log"
  }
}


# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    ${gatk_path} --java-options "-Dsamjdk.compression_level=5 ${java_opt}" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${output_bam_basename}.bam \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER "unsorted" \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 2000000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}
