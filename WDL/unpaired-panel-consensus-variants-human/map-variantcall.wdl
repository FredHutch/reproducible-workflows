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

  # Map reads to reference
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
}

# TASK DEFINITIONS

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
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
    set -o pipefail
    set -e

    java -Dsamjdk.compression_level=5 -Xms3000m -jar /usr/gitc/picard.jar \
      SamToFastq \
			INPUT=${input_bam} \
			FASTQ=/dev/stdout \
			INTERLEAVE=true \
			NON_PF=true | /usr/gitc/bwa mem -K 100000000 -p -v 3 -t 16 -Y ${ref_fasta} /dev/stdin -  2> >(tee ${base_file_name}.bwa.stderr.log >&2) |
		samtools view -1 - > ${base_file_name}.bam

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
