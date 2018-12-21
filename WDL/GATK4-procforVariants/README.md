# [gatk4-data-processing](https://github.com/gatk-workflows/gatk4-data-processing)
## GATK Notes:
### Purpose :
Workflows for processing high-throughput sequencing data for variant discovery with GATK4 and related tools.

### processing-for-variant-discovery-gatk4 :
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the GATK Best Practices
(June 2016).  

#### Requirements/expectations
- Pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
  - filenames all have the same suffix: ".unmapped.bam"
  - files must pass validation by ValidateSamFile
  - reads are provided in query-sorted order
  - all reads must have an RG tag

#### Outputs
- A clean BAM file and its index, suitable for variant discovery analyses.

### Software version requirements :
- GATK 4 or later
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support
 - Successfully tested on v32
 - Does not work on versions < v23 due to output syntax

## Fred Hutch Adaptation
**Workflow:** wdl,  GATK Best Practices [pre-processing for variant calling.](https://github.com/gatk-workflows/gatk4-data-processing)

**Inputs/Parameters:** `sample.inputs.json` defines hg38 as the reference genome, and indicates the location in AWS S3 where the Broad Resource Bundle reference files for this genome are

**Batch File:** `exampleBatch.txt`, indicating the locations of the unmapped bams (with filenames of `sampleName.unmapped.bam`), the location of the bed file associated with each dataset in S3 in tab separated format with a single header row.  

**Output:** A clean BAM file ready for downstream variant calling and it's index written back to S3.
