## Consensus variant calling workflow for human panel-based DNA sequencing.
## Input requirements:
## - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Processing:
## 
##
##
## Output :
## - TBD
## 
## Software version requirements (see recommended dockers in inputs JSON)
## - GATK 4 or later (see gatk docker)
## - Picard (see gotc docker)
## - Samtools (see gotc docker)
##
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

  String annovar_protocols="refGene,knownGene,dbnsfp35c,cosmic70,esp6500siv2_all,exac03,exac03nontcga,1000g2015aug,avsnp150,clinvar_20180603"
  String annovar_operation="g,f,f,f,f,f,f,f,f,f"
  File annovardb
  File table_annovar


scatter (job in batchInfo){
  String sampleName = job.sampleName
  File bamLocation = job.bamLocation
  File bedLocation = job.bedLocation

# Get the basename, i.e. strip the filepath and the extension
  String bam_basename = basename(bamLocation, ".unmapped.bam")
  String base_file_name = sampleName + "." + ref_name



#  # convert unmapped bam to fastq
#   call SamToFastq {
#     input:
#       input_bam = bamLocation,
#       base_file_name = base_file_name
#   }

#     #  Map reads to reference
#   call BwaMem {
#     input:
#       input_fastq = SamToFastq.output_fastq,
#       base_file_name = base_file_name,
#       ref_fasta = ref_fasta,
#       ref_fasta_index = ref_fasta_index,
#       ref_dict = ref_dict,
#       ref_alt = ref_alt,
#       ref_amb = ref_amb,
#       ref_ann = ref_ann,
#       ref_bwt = ref_bwt,
#       ref_pac = ref_pac,
#       ref_sa = ref_sa
#   }

#   # Merge original uBAM and BWA-aligned BAM
#   call MergeBamAlignment {
#     input:
#       unmapped_bam = bamLocation,
#       aligned_bam = BwaMem.output_bam,
#       base_file_name = base_file_name,
#       ref_fasta = ref_fasta,
#       ref_fasta_index = ref_fasta_index,
#       ref_dict = ref_dict
#   }




 # convert unmapped bam to fastq, Map reads to reference
  call SamToFastqAndBwaMem {
    input:
      input_bam = bamLocation,
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
        bed_file = bedLocation,
        base_file_name = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf
    }


    # Generate bcftools vcf
    call bcftoolsMpileup {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        bed_file = bedLocation,
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
        input_SAM_vcf = bcftoolsMpileup.output_vcf,
        ref_name = ref_name,
        base_file_name = base_file_name,
        annovardb = annovardb,
        annovar_protocols = annovar_protocols,
        annovar_operation = annovar_operation,
        table_annovar = table_annovar
    }

  # End scatter 
  }
  # Outputs that will be retained when execution is complete
  output {
    Array[File] bqsr_report = BaseRecalibrator.recalibration_report
    Array[File] analysis_ready_bam = ApplyBQSR.recalibrated_bam 
    Array[File] analysis_ready_bai = ApplyBQSR.recalibrated_bai
    Array[File] analysis_ready_bam_md5 = ApplyBQSR.recalibrated_bam_md5
    Array[File] GATK_vcf = HaplotypeCaller.output_vcf
    Array[File] SAM_vcf = bcftoolsMpileup.output_vcf
    Array[File] GATK_annotated = annovarConsensus.output_GATK_annotated
    Array[File] SAM_annotated = annovarConsensus.output_SAM_annotated
  }
# End workflow
}

# TASK DEFINITIONS


# # Read unmapped BAM, convert to FASTQ
# task SamToFastq {
#   File input_bam
#   String base_file_name

#   command {
#     set -e
#     set -o pipefail

#     java -Dsamjdk.compression_level=5 -Xms3000m -jar /opt/conda/share/picard-2.3.0-0/picard.jar \
#       SamToFastq \
# 			INPUT=${input_bam} \
# 			FASTQ=${base_file_name}.fastq \
# 			INTERLEAVE=true \
# 			NON_PF=true 
#   }
#   runtime {
#     docker: "biocontainers/picard:v2.3.0_cv3"
#     memory: "14 GB"
#     cpu: "1"
#   }
#   output {
#     File output_fastq = "${base_file_name}.fastq"
#   }
# }

# # align to genome
# task BwaMem {
#   File input_fastq
#   String base_file_name
#   File ref_fasta
#   File ref_fasta_index
#   File ref_dict
#   File? ref_alt
#   File ref_amb
#   File ref_ann
#   File ref_bwt
#   File ref_pac
#   File ref_sa

#   command {
#     set -e
#     set -o pipefail

#     ls /cromwell_root/fh-ctr-public-reference-data/genome_data/human/hg38/
    
#     bwa mem \
#       -p -v 3 -t 16 -M \
#       ${ref_fasta} ${input_fastq} | samtools view -1 -bS > ${base_file_name}.aligned.bam 

#   }
#   runtime {
#     docker: "biocontainers/bwa:v0.7.15_cv3"
#     memory: "14 GB"
#     cpu: "16"
#   }
#   output {
#     File output_bam = "${base_file_name}.aligned.bam"
#     File output_bwa_error = "${base_file_name}.bwa.stderr.log"
#   }
# }

# Read unmapped BAM, convert to FASTQ, align to genome
task SamToFastqAndBwaMem {
  File input_bam
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
    set -e
    set -o pipefail

    ls /cromwell_root/fh-ctr-public-reference-data/genome_data/human/hg38/
    
    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=${base_file_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true 

    /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${base_file_name}.fastq | samtools view -1 -bS > ${base_file_name}.aligned.bam 

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14G"
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
    memory: "8G"
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
      --OUTPUT sorted.bam \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false

    /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" \
      SetNmAndUqTags \
      --INPUT sorted.bam \
      --OUTPUT ${base_file_name}.sorted.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}

  }

  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "8G"
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
    memory: "8G"
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
      -bqsr ${recalibration_report} \
      -I ${input_bam} \
      -O ${base_file_name}.recal.bam \
      -R ${ref_fasta} \
      -L ${baserecal_bed_file} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 
  }
  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "8G"
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
    memory: "14G"
    cpu: "1"
  }

  output {
    File output_vcf = "${base_file_name}.GATK.vcf"
    File output_vcf_index = "${base_file_name}.GATK.vcf.tbi"
  }
}



# bcftools Mpileup variant calling
task bcftoolsMpileup {
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
    docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    memory: "8G"
    cpu: "1"
  }

  output {
    File output_vcf = "${base_file_name}.SAM.vcf"
  }
}


# annotate with annovar
task annovarConsensus {
  File input_GATK_vcf
  File input_SAM_vcf
  String base_file_name
  String ref_name
  File annovardb
  String annovar_operation
  String annovar_protocols
  File table_annovar

  command {
  set -e
  
  tar -C annovar/ -xvf ${annovardb}
  
  perl ${table_annovar} ${input_GATK_vcf} annovar/ \
    -buildver ${ref_name} \
    -outfile ${base_file_name}.GATK \
    -remove \
    -protocol ${annovar_protocols} \
    -operation ${annovar_operation} \
    -nastring . -vcfinput -csvout

  perl ${table_annovar} ${input_SAM_vcf} annovar/ \
    -buildver ${ref_name} \
    -outfile ${base_file_name}.SAM \
    -remove \
    -protocol ${annovar_protocols} \
    -operation ${annovar_operation} \
    -nastring . -vcfinput -csvout

  rmdir -r annovar 
  }

  runtime {
    docker: "perl:5.28.0"
    memory: "2G"
    cpu: "1"
  }

  output {
    File output_GATK_annotated = "${base_file_name}.GATK.{ref_name}_multianno.csv"
    File output_SAM_annotated = "${base_file_name}.SAM.{ref_name}_multianno.csv"
  }
}