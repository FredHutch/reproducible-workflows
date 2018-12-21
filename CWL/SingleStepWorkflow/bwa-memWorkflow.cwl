cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
inputs:

  # BWA parameters
  fq1: File
  fq2: File
  genome_index: File
  mem: boolean
  mark_shorter_split_hits: boolean
  read_group_header_line: string
  thread: int?


outputs:

  bwa_result:
    type: File
    outputSource: bwa/tmpsam


steps:
  bwa:
    run: bwa-pe.cwl
    in:
      fq1: fq1
      fq2: fq2
      genome_index: genome_index
      mem: mem
      mark_shorter_split_hits: mark_shorter_split_hits
      read_group_header_line: read_group_header_line
      process: thread
    out: [tmpsam]
