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

  # convert unmapped bam to fastq, Map reads to reference
  call SamToFastqAndBwaMem {
    input:
      input_bam = bamLocation,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  # Merge original uBAM and BWA-aligned BAM
  call MergeBamAlignment {
    input:
      unmapped_bam = bamLocation,
      aligned_bam = SamToFastqAndBwaMem.output_bam,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  # Sort aggregated BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MergeBamAlignment.output_bam,
      base_file_name = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  # Generate the recalibration model by interval
  call BaseRecalibrator {
    input:
      input_bam = SortAndFixTags.output_bam,
      input_bam_index = SortAndFixTags.output_bam_index,
      base_file_name = base_file_name,
      baserecal_bed_file = bedLocation,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
    }


    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        base_file_name = base_file_name,
        recalibration_report = BaseRecalibrator.recalibration_report,
        baserecal_bed_file = bedLocation,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }

    # Generate haplotype caller vcf
    call HaplotypeCaller {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        interval_list = bedLocation,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf
    }


    # Generate samtools vcf
    call samtoolsMpileup {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        interval_list = bedLocation,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf
    }


    # Annotate both sets of variants
    call annovarConsensus {
      input:
        input_GATK_vcf = HaplotypeCaller.output_vcf,
        input_SAM_vcf = samtoolsMpileup.output_vcf,
        ref_name = ref_name
    }



  # Outputs that will be retained when execution is complete
  output {
    #File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = BaseRecalibrator.recalibration_report
    File analysis_ready_bam = ApplyBQSR.recalibrated_bam 
    File analysis_ready_bai = ApplyBQSR.recalibrated_bai
    File analysis_ready_bam_md5 = ApplyBQSR.recalibrated_bam_md5
    File GATK_vcf = HaplotypeCaller.output_vcf
    File Sam_vcf = samtoolsMpileup.output_vcf
    File GATK_annotated = annovarConsensus.output_GATK_annotated
    File SAM_annotated = annovarConsensus.output_SSAM_annotated
    # Add in data provenance outputs
  }

  # End scatter 
  }
# End workflow
}

# TASK DEFINITIONS

# Read unmapped BAM, convert to FASTQ, align to genome
task SamToFastqAndBwaMem {
  File input_bam
  String base_file_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command {
    set -e

    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=${base_file_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true 

    /usr/gitc/bwa index ${ref_fasta}

    /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${base_file_name}.fastq | samtools view \
      -@16 -b -1 > ${base_file_name}.bam \
      2> ${base_file_name}.bwa.stderr.log

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14 GB"
    cpu: "16"
  }
  output {
    File output_fastq = "${base_file_name}.fastq"
    File output_bam = "${base_file_name}.aligned.bam"
    File output_bwa_error = "${base_file_name}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  String base_file_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command {
    # set the bash variable needed for the command-line
    /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${aligned_bam} \
      --UNMAPPED_BAM ${unmapped_bam} \
      --OUTPUT ${base_file_name}.merged.bam \
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
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true 
  }
  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "3500 MB"
    cpu: "1"
  }
  output {
    File output_bam = "${base_file_name}.merged.bam"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String base_file_name
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    set -o pipefail

    /gatk/gatk --java-options "-Dsamjdk.compression_level=$5 -Xms4000m" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /cromwell_root/atempfile \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false

   /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" \
      SetNmAndUqTags \
      --INPUT /cromwell_root/atempfile \
      --OUTPUT ${base_file_name}.sorted.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  rm /cromwell_root/atempfile
  }

  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "5000 MB"
    cpu: "1"
  }
  output {
    File output_bam = "${base_file_name}.sorted.bam"
    File output_bam_index = "${base_file_name}.sorted.bai"
    File output_bam_md5 = "${base_file_name}.sorted.bam.md5"
  }
}



# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File baserecal_bed_file  
  File input_bam_index
  String base_file_name
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    /gatk/gatk --java-options "-Xms4000m" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${base_file_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
       -L ${baserecal_bed_file}
  }
  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "6 GB"
    cpu: "1"
  }
  output {
    File recalibration_report = "${base_file_name}.recal_data.csv"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String base_file_name
  File recalibration_report
  File baserecal_bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {

    /gatk/gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -OBI true \
      -OBM true \
      -O ${base_file_name}.recal.bam \
      -L ${baserecal_bed_file} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --use-original-qualities
  }
  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "3500 MB"
    cpu: "1"
  }
  output {
    File recalibrated_bam = "${base_file_name}.recal.bam"
    File recalibrated_bai = "${base_file_name}.recal.bai"
    File recalibrated_bam_md5 = "${base_file_name}.recal.bam.md5"
  }
}


# HaplotypeCaller per-sample
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  String base_file_name
  File bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf

  command {
    set -e
  
    /gatk/gatk --java-options "-Xmx4g" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --dbsnp ${dbSNP_vcf}\
      -L ${bed_file} \
      -O ${base_file_name}.GATK.vcf 
    }

  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "14 GB"
    cpu: "1"
  }

  output {
    File output_vcf = "${base_file_name}.GATK.vcf"
    File output_vcf_index = "${base_file_name}.GATK.vcf.tbi"
  }
}



# Samtools Mpileup variant calling
task samtoolsMpileup {
  File input_bam
  File input_bam_index
  String base_file_name
  File bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf

  command {
    set -e
  
    samtools mpileup -FORMAT AD
    bcftools mpileup \
      -T ${bed_file} \
      -Ov \
      --annotate "FORMAT/AD,FORMAT/DP" \
      --fasta-ref ${ref_fasta} \
      ${input_bam} | \
    bcftools call -Ov -mv \
          -T ${bed_file} \
          -o ${base_file_name}.SAM.vcf

    }

  runtime {
    docker: docker_image
    memory: mem_size
    cpu: 
  }

  output {
    File output_vcf = "${base_file_name}.SAM.vcf"
  }
}


# annotate with annovar
task annovarConsensus {
  File input_GATK_vcf
  File input_SAM_vcf
  String ref_name

  command {
  set -e
  
  ANOVAR_PROTOCOLS="refGene,esp6500siv2_all,exac03,exac03nontcga,cosmic70,1000g2015aug_all, 1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,snp138,dbnsfp30a"

  perl $ANNOVAR/table_annovar.pl QC/${1}.raw.snps.indels.vcf $ANNOVARDB \
    -buildver ${ref_name} \
    -out annovar/${1} \
    -remove \
    -protocol ${ANOVAR_PROTOCOLS} \
    -operation g,f,f,f,f,f,f,f,f,f,f \
    -nastring . -vcfinput

  perl $ANNOVAR/table_annovar.pl QC/${1}.raw.snps.indels.vcf $ANNOVARDB \
    -buildver ${ref_name} \
    -out annovar/${1} \
    -remove \
    -protocol ${ANOVAR_PROTOCOLS} \
    -operation g,f,f,f,f,f,f,f,f,f,f \
    -nastring . -vcfinput
  }

  runtime {
    docker: "athanadd/annovar"
    memory: mem_size
    cpu: "1"
  }

  output {
    File output_GATK_annotated = "${base_file_name}.GATK.annotated.txt"
    File output_SAM_annotated = "${base_file_name}.SAM.annotated.txt"
  }
}