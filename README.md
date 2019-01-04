# Reproducible Workflows for Scientific Data Analysis

The goal of this repository is to collect a set of workflows that can be
used for reproducible scientific data analysis. In general these workflows
are designed to be run using Cromwell, although in principal this repo
could also be used to keep track of workflows for any other language or
system.


### Workflow Folder Organization

Each workflow is stored in its own folder, which contains:
  * The workflow file (e.g. .wdl).
  * An example inputs/parameters file (e.g. .json).
  * An example batch file (e.g. .csv).
  * A readme (in Markdown) explaining what the workflow does, tool descriptions and guidance for inputs/dependencies/outputs and how to run the workflow on the sample data.

### Workflow Design Principals

In developing reproducible workflows the approach one uses to defining the tools used in a workflow, any reference/dependencies required, the particular inputs desired and the data upon which to operate in a given batch can impact their ability to be reused.  Here we describe some guidelines for approaching workflow design.  

In general, we subscribe to the philosophy of data organization that a workflow describes a series of tools/processes that are relatively general.  This workflow can be passed a group of inputs/parameters that modify how the workflow functions where needed (such as the application of a workflow for a research project).  The input data sets and any data set-specific parameters that the researcher intends to apply the same workflow and input/parameters to is described in a batch file.  These three files together along with any reference data and raw input data, should completely describe a process that transforms the input data to the resulting output files for it to be reproducible.  

#### Cromwell/WDL/AWS workflows
1.  Workflow Description file:
  - in WDL, a list of tools to be run in a sequence, likely several, otherwise using a workflow manager is not the right approach.  
  - This file describes the process that is desired to occur every time the workflow is run.
2.  Inputs/Parameters file:
  - in json, a workflow-specific list of inputs and parameters that are intended to be set for every group of workflow executions.
  - Examples of what this input may include would be which genome to map to, reference data files to use, etc.
3.  Batch file:
  - in csv, a batch-specific list of the raw input data sets intended to be processed using the same set of inputs/parameters for the same workflow.  
  - This file is a list of data locations and any other sample/job-specific information the workflow needs.  Ideally this would be relatively minimal so that the consistency of the analysis between input data sets are as similar as possible to leverage the strengths of a reproducible workflow.  


### Workflow Development Status

A validated tool is one that can be successfully executed using the WDL
and inputs JSON, and is determined by the developer to be ready for use.  These must also include a sample data set in a publicly accessible location (e.g., AWS S3).  

#### Validated Workflows

  * WDL / align-proteins-diamond: Aligns a set of proteins using DIAMOND
  * WDL / denovo-gene-level-metagenomics-cags: Assemble a set of metagenomes, quantify genes, and make CAGs

#### Under Development Workflows

  * Coming soon...


### Useful Workflow Development tools and resources

https://cromwell.readthedocs.io/en/develop/WOMtool/

Web based WDL viewer
http://pb.opensource.epam.com/
