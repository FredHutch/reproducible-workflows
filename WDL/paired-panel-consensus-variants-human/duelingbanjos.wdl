## Consensus variant calling workflow for human panel-based DNA sequencing for paired samples.
## Input requirements:
## - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference bam and sample bam must be provided
## - bed file for the panel used to generate both datasets must be provided (and must be the same one)
## 
## Processing:
## 
##
##
## Output :
## - recalibrated bam and it's index for the sample and the reference
## - 
## - GATK/mutect2 vcf
## - strelka2 vcf
## - Annovar annotated vcfs and tabular file for each variant caller
## 

## Software version requirements (see recommended dockers in inputs JSON)
## - GATK 4 or later (see gatk docker)
## - Picard (see gotc docker)
## - Samtools (see gotc docker)
##
workflow Panel_BWA_GATK4_Strelka_Var_Annotate_Paired {
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

  File annovarTAR
  String annovar_protocols
  String annovar_operation


scatter (job in batchInfo){
  String sampleName = job.sampleName
  File refBamLocation = job.refBamLocation
  File sampleBamLocation = job.sampleBamLocation
  File bedLocation = job.bedLocation

# Get the basename, i.e. strip the filepath and the extension
  String bam_basename = basename(sampleBamLocation, ".unmapped.bam")
  String ref_basename = basename(refBamLocation, ".unmapped.bam")
  String base_file_name = sampleName + "." + ref_name
  String ref_file_name = ref_basename + "." + ref_name


# convert unmapped bam to fastq
  call SamToFastq {
    input:
      input_bam = sampleBamLocation,
      ref_bam = refBamLocation,
      base_file_name = base_file_name,
      ref_file_name = ref_file_name
  }

#  Map reads to reference
  call BwaMem {
    input:
      input_fastq = SamToFastq.output_fastq,
      ref_fastq = SamToFastq.output_ref_fastq,
      base_file_name = base_file_name,
      ref_file_name = ref_file_name,
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
      unmapped_bam = sampleBamLocation,
      unmapped_ref_bam = refBamLocation,
      aligned_bam = BwaMem.output_bam,
      aligned_ref_bam = BwaMem.output_ref_bam,
      base_file_name = base_file_name,
      ref_file_name = ref_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }

  # Generate the recalibration model by interval
  call BaseRecalibrator {
    input:
      input_bam = MergeBamAlignment.output_bam,
      input_bam_index = MergeBamAlignment.output_bai,
      base_file_name = base_file_name,
      input_ref_bam = MergeBamAlignment.output_ref_bam,
      input_ref_bam_index = MergeBamAlignment.output_ref_bai,
      ref_file_name = ref_file_name,
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
        input_bam = MergeBamAlignment.output_bam,
        input_bam_index = MergeBamAlignment.output_bai,
        base_file_name = base_file_name,
        input_ref_bam = MergeBamAlignment.output_ref_bam,
        input_ref_bam_index = MergeBamAlignment.output_ref_bai,
        ref_file_name = ref_file_name,
        recalibration_report = BaseRecalibrator.recalibration_report,
        recalibration_report_ref = BaseRecalibrator.recalibration_report_ref,
        baserecal_bed_file = bedLocation,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }

    # Generate Mutect2 vcf
    call Mutect2 {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        base_file_name = base_file_name,
        input_ref_bam = ApplyBQSR.recalibrated_ref_bam,
        input_ref_bam_index = ApplyBQSR.recalibrated_ref_bai,
        ref_file_name = ref_file_name,
        bed_file = bedLocation,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf
    }


    # Generate strelka2 vcf
    call Strelka2 {
      input:
        input_bam = ApplyBQSR.recalibrated_bam,
        input_bam_index = ApplyBQSR.recalibrated_bai,
        base_file_name = base_file_name,
        input_ref_bam = ApplyBQSR.recalibrated_ref_bam,
        input_ref_bam_index = ApplyBQSR.recalibrated_ref_bai,
        ref_file_name = ref_file_name,
        bed_file = bedLocation,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        dbSNP_vcf = dbSNP_vcf
    }

    # concatenate the VCF's from Strelka
    call ConcatVCFs {
      input:
      first_vcf = Strelka2.output_indels_vcf,
      second_vcf = Strelka2.output_snvs_vcf
    }

    # Annotate both sets of variants
    call annovarConsensus {
      input:
        input_mutect_vcf = Mutect2.output_vcf,
        input_strelka_vcf = ConcatVCFs.merged_vcf,
        ref_name = ref_name,
        base_file_name = base_file_name,
        annovarTAR = annovarTAR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols
    }

  # End scatter 
  }
  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysis_ready_bam = ApplyBQSR.recalibrated_bam 
    Array[File] analysis_ready_bai = ApplyBQSR.recalibrated_bai
    Array[File] analysis_ready_ref_bam = ApplyBQSR.recalibrated_ref_bam 
    Array[File] analysis_ready_ref_bai = ApplyBQSR.recalibrated_ref_bai
    Array[File] output_re_bam = Mutect2.output_re_bam
    Array[File] mutect_vcf = Mutect2.output_vcf
    Array[File] strelka_vcf = ConcatVCFs.merged_vcf
    Array[File] mutect_annotated_vcf = annovarConsensus.output_mutect_annotated_vcf
    Array[File] mutect_annotated = annovarConsensus.output_mutect_annotated_table
    Array[File] strelka_annotated_vcf = annovarConsensus.output_strelka_annotated_vcf
    Array[File] strelka_annotated = annovarConsensus.output_strelka_annotated_table
  }
# End workflow
}

# TASK DEFINITIONS


# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  File input_bam
  File ref_bam
  String base_file_name
  String ref_file_name

  command {
    set -eo pipefail

    find /cromwell_root/ -type f
    
    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=${base_file_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true 

    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${ref_bam} \
			FASTQ=${ref_file_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true 
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14 GB"
    cpu: "4"
  }
  output {
    File output_fastq = "${base_file_name}.fastq"
    File output_ref_fastq = "${ref_file_name}.fastq"
  }
}

# align to genome
task BwaMem {
  File input_fastq
  String base_file_name
  File ref_fastq
  String ref_file_name
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
    set -eo pipefail
    
    /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${input_fastq} | samtools view -1bS > ${base_file_name}.aligned.bam 

    /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${ref_fastq} | samtools view -1bS > ${ref_file_name}.aligned.bam 

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    memory: "14 GB"
    cpu: "16"
  }
  output {
    File output_bam = "${base_file_name}.aligned.bam"
    File output_ref_bam = "${ref_file_name}.aligned.bam"
  }
}


# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  String base_file_name
  File unmapped_ref_bam
  File aligned_ref_bam
  String ref_file_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command {
    set -eo pipefail

    /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      MergeBamAlignment \
     --ALIGNED_BAM ${aligned_bam} \
     --UNMAPPED_BAM ${unmapped_bam} \
     --OUTPUT ${base_file_name}.merged.bam \
     --REFERENCE_SEQUENCE ${ref_fasta} \
     --PAIRED_RUN true \
     --SORT_ORDER coordinate \
     --CREATE_INDEX true \
     --CLIP_ADAPTERS true \
     --MAX_RECORDS_IN_RAM 2000000 \
     --MAX_INSERTIONS_OR_DELETIONS -1 \
     --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
     --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
     --ALIGNER_PROPER_PAIR_FLAGS true

    /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      MergeBamAlignment \
     --ALIGNED_BAM ${aligned_ref_bam} \
     --UNMAPPED_BAM ${unmapped_ref_bam} \
     --OUTPUT ${ref_file_name}.merged.bam \
     --REFERENCE_SEQUENCE ${ref_fasta} \
     --PAIRED_RUN true \
     --SORT_ORDER coordinate \
     --CREATE_INDEX true \
     --CLIP_ADAPTERS true \
     --MAX_RECORDS_IN_RAM 2000000 \
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
    File output_bai = "${base_file_name}.merged.bai"
    File output_ref_bam = "${ref_file_name}.merged.bam"
    File output_ref_bai = "${ref_file_name}.merged.bai"
  }
}



# Generate Base Quality Score Recalibration (BQSR) models
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String base_file_name
  File input_ref_bam
  File input_ref_bam_index
  String ref_file_name
  File baserecal_bed_file  
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    set -eo pipefail

    /gatk/gatk --java-options "-Xms4g" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${base_file_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
       -L ${baserecal_bed_file}    
    
    /gatk/gatk --java-options "-Xms4g" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_ref_bam} \
      --use-original-qualities \
      -O ${ref_file_name}.recal_data.csv \
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
    File recalibration_report_ref = "${ref_file_name}.recal_data.csv"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String base_file_name
  File input_ref_bam
  File input_ref_bam_index
  String ref_file_name
  File recalibration_report
  File recalibration_report_ref
  File baserecal_bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    set -eo pipefail

    /gatk/gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -bqsr ${recalibration_report} \
      -I ${input_bam} \
      -O ${base_file_name}.recal.bam \
      -R ${ref_fasta} \

    /gatk/gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -bqsr ${recalibration_report_ref} \
      -I ${input_ref_bam} \
      -O ${ref_file_name}.recal.bam \
      -R ${ref_fasta} \
  }
  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "8G"
    cpu: "1"
  }
  output {
    File recalibrated_bam = "${base_file_name}.recal.bam"
    File recalibrated_bai = "${base_file_name}.recal.bai"
    File recalibrated_ref_bam = "${ref_file_name}.recal.bam"
    File recalibrated_ref_bai = "${ref_file_name}.recal.bai"
  }
}


# Mutect 2 calling
task Mutect2 {
  File input_bam
  File input_bam_index
  String base_file_name
  File input_ref_bam
  File input_ref_bam_index
  String ref_file_name
  File bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf

  command {
    set -eo pipefail
  
    /gatk/gatk --java-options "-Xmx2g" \
      Mutect2 \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -I ${input_ref_bam} \
      -normal ${ref_file_name} \
      -tumor ${base_file_name} \
      -O ${base_file_name}.mutect2.vcf 
      -bamout ${base_file_name}_${ref_file_name}.reassembled.bam
    }

  runtime {
    docker: "broadinstitute/gatk:4.0.4.0"
    memory: "14G"
    cpu: "4"
  }

  output {
    File output_vcf = "${base_file_name}.mutect2.vcf"
    File output_vcf_index = "${base_file_name}.mutect2.vcf.tbi"
    File output_re_bam = "${base_file_name}_${ref_file_name}.reassembled.bam"
  }
}



# strelka2 variant calling FINISH TRANSITION FROM SAMTOOLS
task Strelka2 {
  File input_bam
  File input_bam_index
  String base_file_name
  File input_ref_bam
  File input_ref_bam_index
  String ref_file_name
  File bed_file
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File dbSNP_vcf

  command {
    set -eo pipefail
  

  configureStrelkaSomaticWorkflow.py \
    --normalBam=${input_ref_bam} \
    --tumorBam=${input_bam}\
    --referenceFasta=${ref_fasta} \
    --targeted \
    --runDir=strelkatemp
  
  strelkatemp/runWorkflow.py
  gunzip strelkatemp/results/variants/somatic.snvs.vcf.gz ${base_file_name}.somatic.snvs.vcf 
  gunzip strelkatemp/results/variants/somatic.indels.vcf.gz ${base_file_name}.somatic.indels.vcf 

    }

  runtime {
    docker: "quay.io/biocontainers/strelka:2.9.10--0"
    memory: "14G"
    cpu: "4"
  }

  output {
    File output_snvs_vcf = "${base_file_name}.somatic.snvs.vcf"
    File output_indels_vcf = "${base_file_name}.somatic.indels.vcf"
}
##--callRegions=${bed_file}  because it needs a bzip'd tabix indexed bed file bgzip -c file.vcf | tabix -p vcf bedfile.vcf.gz to get vcf.gz and vcf.gz.tbi
## https://github.com/BioContainers/containers/blob/master/strelka/2.9.7/Dockerfile
}


# Concatenate STrelka vcf's
task ConcatVCFs {
  File first_vcf
  File second_vcf

  command {
    set -eo pipefail

    bcftools sort -O v -o first.vcf ${first_vcf}
    bcttools sort -O v -o second.vcf ${second_vcf}
    bcftools concat -a -O v -o strelka.merged.vcf first.vcf second.vcf 
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    memory: "2G"
    cpu: "1"
  }
  output {
    File merged_vcf = "strelka.merged.vcf"
  }
}

# annotate with annovar
task annovarConsensus {
  File input_mutect_vcf
  File input_strelka_vcf
  String base_file_name
  String ref_name
  File annovarTAR
  String annovar_protocols
  String annovar_operation

  command {
   set -eo pipefail
  
   tar -xzvf ${annovarTAR}
  
    perl annovar/table_annovar.pl ${input_mutect_vcf} annovar/humandb/ \
      -buildver ${ref_name} \
      -outfile ${base_file_name}.GATK \
      -remove \
     -protocol ${annovar_protocols} \
     -operation ${annovar_operation} \
      -nastring . -vcfinput

   perl annovar/table_annovar.pl ${input_strelka_vcf} annovar/humandb/ \
      -buildver ${ref_name} \
      -outfile ${base_file_name}.strelka \
      -remove \
      -protocol ${annovar_protocols} \
      -operation ${annovar_operation} \
      -nastring . -vcfinput
  }

  runtime {
    docker: "perl:5.28.0"
    memory: "4G"
    cpu: "1"
  }

  output {
    File output_mutect_annotated_vcf = "${base_file_name}.mutect.${ref_name}_multianno.vcf"
    File output_mutect_annotated_table = "${base_file_name}.mutect.${ref_name}_multianno.txt"
    File output_strelka_annotated_vcf = "${base_file_name}.strelka.${ref_name}_multianno.vcf"
    File output_strelka_annotated_table = "${base_file_name}.strelka.${ref_name}_multianno.txt"
  }
}