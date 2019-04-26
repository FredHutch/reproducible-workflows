#!/usr/bin/env nextflow

// Set default values for parameters
// These can be changed when invoking the script by setting (e.g.) --ref_name hg37
params.ref_name = "hg38"
ref_name = params.ref_name

// REFERENCE FILES
params.ref_fasta = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta"
params.ref_fasta_index = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.fai"
params.ref_dict = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.dict"
params.ref_pac = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.pac"
params.ref_sa = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.sa"
params.ref_amb = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.amb"
params.ref_ann = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.ann"
params.ref_alt = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.alt"
params.ref_bwt = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta.bwt"
params.dbSNP_vcf_index = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
params.dbSNP_vcf = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"

// ANNOVAR TOOL AND ASSOCIATED DATABASE FILES
params.annovarTAR = "s3://fh-ctr-public-reference-data/tool_specific_data/annovar/2018Apr16-ANNOVAR-2019Jan15-hg38humandb.tar.gz"
params.annovar_protocols = "refGene,knownGene,dbnsfp35c,cosmic70,esp6500siv2_all,exac03,exac03nontcga,avsnp150,clinvar_20180603"
params.annovar_operation = "g,f,f,f,f,f,f,f,f"

// Specify files as comma-delimited string
params.known_indels_sites_indices = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi,s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
known_indels_sites_indices_ch = Channel.from(params.known_indels_sites_indices)
                                       .map { it -> it.split(",") }
                                       .flatMap()
                                       .map { it -> file(it) }

params.known_indels_sites_VCFs = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
known_indels_sites_VCFs_ch = Channel.from(params.known_indels_sites_VCFs)
                                    .map { it -> it.split(",") }
                                    .flatMap()
                                    .map { it -> file(it) }

// Two identical channels, each associating samples with BAM files
Channel.from(file(params.batchfile))
       .splitCsv(header: true, sep: "\t")
       .map { job ->
       [job.sampleName, file(job.bamLocation)]}
       .into{ bam_for_SamToFastq_ch; bam_for_MergeBamAlignment_ch }

// Channel linking each sample to its BED file
Channel.from(file(params.batchfile))
       .splitCsv(header: true, sep: "\t")
       .map { job ->
       [job.sampleName, file(job.bedLocation)]}
       .into{ baserecal_bed_BaseRecalibrator_ch; baserecal_bed_ApplyBQSR_ch; baserecal_bed_HaplotypeCaller_ch; baserecal_bed_Mpileup_ch }


// First task, which consumes each sample
process SamToFastq {
    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    cpus 4
    memory "14 GB"
    // errorStrategy "retry"

    input:
    set val(sampleName), file(input_bam) from bam_for_SamToFastq_ch
    val ref_name from params.ref_name

    output:
    set val("${sampleName}"), file("${sampleName}.${ref_name}.fastq") into fastq_ch

    """
    set -eo pipefail
    
    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=${sampleName}.${ref_name}.fastq \
			INTERLEAVE=true \
			NON_PF=true 
    """
}

// Map reads to reference
process BwaMem {
    container "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    cpus 16
    memory "14 GB"

    input:
    set val(sampleName), file(input_fastq) from fastq_ch
    val ref_name from params.ref_name
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)
    file ref_dict from file(params.ref_dict)
    file ref_alt from file(params.ref_alt)
    file ref_amb from file(params.ref_amb)
    file ref_ann from file(params.ref_ann)
    file ref_bwt from file(params.ref_bwt)
    file ref_pac from file(params.ref_pac)
    file ref_sa from file(params.ref_sa)

    output:
    set val(sampleName), file("${sampleName}.${ref_name}.aligned.bam") into aligned_bam_ch

    """
    set -eo pipefail
    
    /usr/gitc/bwa mem \
      -p -v 3 -t 16 -M \
      ${ref_fasta} ${input_fastq} | samtools view -1bS > ${sampleName}.${ref_name}.aligned.bam 
    """

}

//   Merge original uBAM and BWA-aligned BAM
process MergeBamAlignment {
    container "broadinstitute/gatk:4.0.4.0"
    cpus 1
    memory "8G"

    input:
    set val(sampleName), file(aligned_bam), file(unmapped_bam) from aligned_bam_ch.join(bam_for_MergeBamAlignment_ch)
    val ref_name from params.ref_name
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)
    file ref_dict from file(params.ref_dict)

    output:
    set val("${sampleName}"), file("${sampleName}.${ref_name}.merged.bam"), file("${sampleName}.${ref_name}.merged.bai") into merged_bam_BaseRecalibrator_ch, merged_bam_ApplyBQSR_ch, merged_bam_HaplotypeCaller_ch, merged_bam_Mpileup_ch
        
    """
    set -eo pipefail

    # set the bash variable needed for the command-line
    /gatk/gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      MergeBamAlignment \
     --ALIGNED_BAM ${aligned_bam} \
     --UNMAPPED_BAM ${unmapped_bam} \
     --OUTPUT ${sampleName}.${ref_name}.merged.bam \
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
    """
  }

process BaseRecalibrator {
    container "broadinstitute/gatk:4.0.4.0"
    cpus 1
    memory "8G"

    input:
    set sampleName, file(input_bam), file(input_bam_index), file(baserecal_bed_file) from merged_bam_BaseRecalibrator_ch.join(baserecal_bed_BaseRecalibrator_ch)
    val ref_name from params.ref_name
    file dbSNP_vcf from file(params.dbSNP_vcf)
    file dbSNP_vcf_index from file(params.dbSNP_vcf_index)
    file ref_dict from file(params.ref_dict)
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)
    file known_indels_sites_indices from known_indels_sites_indices_ch.collect()
    file known_indels_sites_VCFs from known_indels_sites_VCFs_ch.collect()

    output:
    set val("${sampleName}"), file("${sampleName}.${ref_name}.recal_data.csv") into recalibration_report_ch

    """
    set -eo pipefail

    /gatk/gatk --java-options "-Xms4g" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${sampleName}.${ref_name}.recal_data.csv \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${known_indels_sites_VCFs.join(" --known-sites ")} \
       -L ${baserecal_bed_file}
    """
}

// Apply Base Quality Score Recalibration (BQSR) model
process ApplyBQSR {
    container "broadinstitute/gatk:4.0.4.0"
    cpus 1
    memory "8G"
    publishDir "${params.output_directory}"

    input:
    set val(sample_name), file(input_bam), file(input_bam_index), file(recalibration_report), file(baserecal_bed_file) from merged_bam_ApplyBQSR_ch.join(recalibration_report_ch).join(baserecal_bed_ApplyBQSR_ch)
    val ref_name from params.ref_name
    file ref_dict from file(params.ref_dict)
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)

    output:
    file "${sample_name}.${ref_name}.recal.bam"
    file "${sample_name}.${ref_name}.recal.bai"

    """
    set -eo pipefail

    /gatk/gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -bqsr ${recalibration_report} \
      -I ${input_bam} \
      -O ${sample_name}.${ref_name}.recal.bam \
      -R ${ref_fasta}

    """
}

// HaplotypeCaller per-sample
process HaplotypeCaller {
    container "broadinstitute/gatk:4.0.4.0"
    cpus 4
    memory "14G"
    publishDir "${params.output_directory}"
    

    input:
    val ref_name from params.ref_name
    set val(sample_name), file(input_bam), file(input_bam_index), file(bed_file) from merged_bam_HaplotypeCaller_ch.join(baserecal_bed_HaplotypeCaller_ch)
    file ref_dict from file(params.ref_dict)
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)
    file dbSNP_vcf from file(params.dbSNP_vcf)
    file dbSNP_vcf_index from file(params.dbSNP_vcf_index)
    
    output:
    set val("${sample_name}"), file("${sample_name}.${ref_name}.GATK.vcf"), file("${sample_name}.${ref_name}.GATK.vcf.tbi") into annovar_ch

    """
    set -eo pipefail

    /gatk/gatk --java-options "-Xmx4g" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${sample_name}.${ref_name}.GATK.vcf 
    """

}


// bcftools Mpileup variant calling
process bcftoolsMpileup {
    container "quay.io/biocontainers/bcftools:1.9--h4da6232_0"
    cpus 4
    memory "14G"
    publishDir "${params.output_directory}"

    input:
    val ref_name from params.ref_name
    set val(sample_name), file(input_bam), file(input_bam_index), file(bed_file) from merged_bam_Mpileup_ch.join(baserecal_bed_Mpileup_ch)
    file ref_dict from file(params.ref_dict)
    file ref_fasta from file(params.ref_fasta)
    file ref_fasta_index from file(params.ref_fasta_index)
    file dbSNP_vcf from file(params.dbSNP_vcf)
    file dbSNP_vcf_index from file(params.dbSNP_vcf_index)
 
    output:
    set val("${sample_name}"), file("${sample_name}.${ref_name}.SAM.vcf") into bcftools_ch


    """
    set -eo pipefail
  
    bcftools mpileup \
      --max-depth 10000 \
      --max-idepth 10000 \
      --annotate "FORMAT/AD,FORMAT/DP" \
      --fasta-ref ${ref_fasta} \
      --ignore-RG \
      --no-BAQ \
      ${input_bam} | bcftools call -Ov -mv \
          -o ${sample_name}.${ref_name}.SAM.vcf

    """
}

// Annotate with annovar
process annovarConsensus {
    container "perl:5.28.0"
    cpus 1
    memory "4G"
    publishDir "${params.output_directory}"

    input:
    set val(sample_name), file(input_GATK_vcf), file(input_GATK_vcf_tbi), file(input_SAM_vcf) from annovar_ch.join(bcftools_ch)
    val ref_name from params.ref_name
    file annovarTAR from file(params.annovarTAR)
    val annovar_protocols from params.annovar_protocols
    val annovar_operation from params.annovar_operation

    output:
    file "${sample_name}.${ref_name}.GATK.${ref_name}_multianno.vcf"
    file "${sample_name}.${ref_name}.GATK.${ref_name}_multianno.txt"
    file "${sample_name}.${ref_name}.SAM.${ref_name}_multianno.vcf"
    file "${sample_name}.${ref_name}.SAM.${ref_name}_multianno.txt"

    """
    set -eo pipefail
  
    tar -xzvf ${annovarTAR}
  
    perl annovar/table_annovar.pl ${input_GATK_vcf} annovar/humandb/ \
      -buildver ${ref_name} \
      -outfile ${sample_name}.${ref_name}.GATK \
      -remove \
      -protocol ${annovar_protocols} \
      -operation ${annovar_operation} \
      -nastring . -vcfinput

   perl annovar/table_annovar.pl ${input_SAM_vcf} annovar/humandb/ \
      -buildver ${ref_name} \
      -outfile ${sample_name}.${ref_name}.SAM \
      -remove \
      -protocol ${annovar_protocols} \
      -operation ${annovar_operation} \
      -nastring . -vcfinput

    """
}