class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: metaphlan2
baseCommand:
  - metaphlan2.py
inputs:
  - id: input_fastq
    type: File
    inputBinding:
      position: 4
    label: Input
    doc: Input file in FASTQ format (gzip supported)
  - id: output_file
    type: string
    label: OutputPath
    doc: Output file of MetaPhlAn2 classifications
outputs:
  - id: output
    doc: Output file created by MetaPhlAn2
    label: Output
    type: File
    outputBinding:
      glob: $(inputs.output_file)
    secondaryFiles: []
label: metaphlan2
arguments:
  - position: 0
    prefix: '--input_type'
    valueFrom: fastq
  - position: 1
    prefix: '--tmp_dir'
    valueFrom: ./
  - position: 2
    prefix: '-o'
    valueFrom: $(inputs.output_file)
  - position: 3
    prefix: '--bowtie2out'
    valueFrom: bowtie2_alignment.txt
requirements:
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 0
  - class: DockerRequirement
    dockerPull: 'quay.io/fhcrc-microbiome/metaphlan:latest'
