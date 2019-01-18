class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: kraken2
baseCommand:
  - kraken2
inputs:
  - id: input_fastq
    type: File
    inputBinding:
      position: 1
    label: Input
    doc: Input file in FASTQ format (gzip supported)
  - id: output_file
    type: string
    label: OutputPath
    doc: Output file of Kraken2 classifications
outputs:
  - id: output
    doc: Output file created by Kraken2
    label: Output
    type: File
    outputBinding:
      glob: $(inputs.output_file)
    secondaryFiles: []
label: kraken2
arguments: [
  "--db",
  "/database",
  "--output",
  "-",
  "--use-mpa-style",
  "--report",
  $(inputs.output_file)
]
requirements:
  - class: ResourceRequirement
    ramMin: 4000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'quay.io/fhcrc-microbiome/kraken2:2.0.7_beta_minikraken2_v1_8GB'
