fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '/*_{1,2}.fastq.gz', flat:true)

process interleave_fastq_pairs {
  publishDir params.output_folder

  container "ubuntu:16.04"
  cpus 1
  memory "1 GB"

  input:
  set pair_name, file(fastq1), file(fastq2) from fastq_pair_ch

  output:
  file "${pair_name}.fastq.gz"

  """
  set -e

  # Some basic checks that the line numbers match
  (( \$(gunzip -c ${fastq1} | wc l) == \$(gunzip -c ${fastq2} | wc l) ))

  # Now interleave the files
  paste <(gunzip -c ${fastq1}) <(gunzip -c ${fastq2}) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "${pair_name}.fastq.gz"
  """
    
}