## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements (see recommended dockers in inputs JSON)
## - GATK 4 or later
## - Picard (see gotc docker)
## - Samtools (see gotc docker)
## - Python 2.7
##
## Cromwell version support
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.

# WORKFLOW DEFINITION
workflow BWA_GATK4_samtools_Variants {
  File batchFile
  Array[Object] batchInfo = read_objects(batchFile)

  String ref_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String bwa_commandline
  Int compression_level

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  String gotc_docker
  String gatk_docker
  String python_docker

  String gotc_path
  String gatk_path

  Int flowcell_small_disk
  Int flowcell_medium_disk
  Int agg_small_disk
  Int agg_medium_disk
  Int agg_large_disk

  Int preemptible_tries
  Int agg_preemptible_tries

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call GetBwaVersion {
    input:
      docker_image = gotc_docker,
      bwa_path = gotc_path,
      preemptible_tries = preemptible_tries
  }

  String BWAVersion = GetBwaVersion.version

scatter (job in batchInfo){
    String sampleName = job.sampleName
    File bamLocation = job.bamLocation
    File bedLocation = job.bedLocation

  # Get the basename, i.e. strip the filepath and the extension
  String bam_basename = basename(bamLocation, ".bam")

  String base_file_name = sampleName + "." + ref_name
  # Get the bed file for this dataset
  # Array[Array[String]] sequenceIntervals = read_tsv(bedLocation)

  # Map reads to reference
  call SamToFastqAndBwaMem {
    input:
      input_bam = bamLocation,
      bwa_commandline = bwa_commandline,
      output_bam_basename = bam_basename + ".unmerged",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      docker_image = gotc_docker,
      bwa_path = gotc_path,
      gotc_path = gotc_path,
      disk_size = flowcell_medium_disk,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  # Merge original uBAM and BWA-aligned BAM
  call MergeBamAlignment {
    input:
      unmapped_bam = bamLocation,
      bwa_commandline = bwa_commandline,
      bwa_version = BWAVersion,
      aligned_bam = SamToFastqAndBwaMem.output_bam,
      output_bam_basename = bam_basename + ".aligned.unsorted",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = flowcell_medium_disk,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }
  # mark duplicates
  # call MarkDuplicates {
  #   input:
  #     input_bams = MergeBamAlignment.output_bam,
  #     output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
  #     metrics_filename = base_file_name + ".duplicate_metrics",
  #     docker_image = gatk_docker,
  #     gatk_path = gatk_path,
  #     disk_size = agg_large_disk,
  #     compression_level = compression_level,
  #     preemptible_tries = agg_preemptible_tries
  # }

  # Sort aggregated BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = agg_large_disk,
      preemptible_tries = 0,
      compression_level = compression_level
  }


# Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        baserecal_bed_file = bedLocation,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }


    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.recalibrated", 
        recalibration_report = BaseRecalibrator.recalibration_report,
        baserecal_bed_file = bedLocation,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        disk_size = agg_small_disk,
        preemptible_tries = agg_preemptible_tries
    }

# Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        interval_list = bedLocation,
        output_filename = base_file_name + ".vcf.gz",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        make_gvcf = false,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }


  # Outputs that will be retained when execution is complete
  output {
    #File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = BaseRecalibrator.recalibration_report
    File analysis_ready_bam = ApplyBQSR.recalibrated_bam 
    File analysis_ready_bai = ApplyBQSR.recalibrated_bai
    File analysis_ready_bam_md5 = ApplyBQSR.recalibrated_bam_md5
    File output_vcf = HaplotypeCaller.output_vcf
    # Add in data provenance outputs
  }
}

# TASK DEFINITIONS

# Get version of BWA
task GetBwaVersion {

  Int preemptible_tries
  String mem_size

  String docker_image
  String bwa_path

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ${bwa_path}bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
  }
  output {
    String version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File input_bam
  String bwa_commandline
  String output_bam_basename
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

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size
  String num_cpu

  String docker_image
  String bwa_path
  String gotc_path
  String java_opt

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

		java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar ${gotc_path}picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true | ${bwa_path}${bwa_commandline} /dev/stdin -  2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) |
		samtools view -1 - > ${output_bam_basename}.bam

  >>>
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
    cpu: num_cpu
    #disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
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
    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
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
    #disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  Int compression_level
  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt_sort
  String java_opt_fix

  command <<<
    set -o pipefail

    ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_sort}" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /cromwell_root/atempfile \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false
   ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt_fix}" \
      SetNmAndUqTags \
      --INPUT /cromwell_root/atempfile \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  rm /cromwell_root/atempfile
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
    #disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# # Mark duplicate reads to avoid counting non-independent observations
# task MarkDuplicates {
#   Array[File] input_bams
#   String output_bam_basename
#   String metrics_filename

#   Int compression_level
#   Int preemptible_tries
#   Int disk_size
#   String mem_size

#   String docker_image
#   String gatk_path
#   String java_opt

#  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
#  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
#  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
#   command {
#     ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
#       MarkDuplicates \
#       --INPUT ${sep=' --INPUT ' input_bams} \
#       --OUTPUT ${output_bam_basename}.bam \
#       --METRICS_FILE ${metrics_filename} \
#       --VALIDATION_STRINGENCY SILENT \
#       --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
#       --ASSUME_SORT_ORDER "queryname" \
#       --CREATE_MD5_FILE true
#   }
#   runtime {
#     preemptible: preemptible_tries
#     docker: docker_image
#     memory: mem_size
#     #disks: "local-disk " + disk_size + " HDD"
#   }
#   output {
#     File output_bam = "${output_bam_basename}.bam"
#     File duplicate_metrics = "${metrics_filename}"
#   }
# }


# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File baserecal_bed_file  
  File input_bam_index
  String recalibration_report_filename
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {

    gunzip ${ref_fasta}
    gunzipped_ref_fasta=$(echo ${ref_fasta} | sed "s/\.gz//")

    ${gatk_path} --java-options "${java_opt}" \
      BaseRecalibrator \
      -R $gunzipped_ref_fasta \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
       -L ${baserecal_bed_file}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
    #disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  File baserecal_bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  Int preemptible_tries
  Int disk_size
  String mem_size

  String docker_image
  String gatk_path
  String java_opt

  command {
    gunzip ${ref_fasta}
    gunzipped_ref_fasta=$(echo ${ref_fasta} | sed "s/\.gz//")

    ${gatk_path} --java-options "${java_opt}" \
      ApplyBQSR \
      -R $gunzipped_ref_fasta \
      -I ${input_bam} \
      -OBI true \
      -OBM true \
      -O ${output_bam_basename}.bam \
      -L ${baserecal_bed_file} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: mem_size
    #disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bai = "${output_bam_basename}.bai"
    File recalibrated_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# HaplotypeCaller per-sample
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String output_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String docker_image
  String gatk_path
  String? java_options
  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

  # Runtime parameters
  Int preemptible_tries
  Int disk_size
  String mem_size


  command <<<
  set -e
  
    ${gatk_path} --java-options "-Xmx${command_mem_gb}G ${java_opt}" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -L ${interval_list} \
      -O ${output_filename} 
  >>>

  runtime {
    docker: docker_image
    preemptible: preemptible_tries
    memory: mem_size
  }

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}