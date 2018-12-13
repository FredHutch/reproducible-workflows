# Reproducible Workflows for Scientific Data Analysis

The goal of this repository is to collect a set of workflows that can be
used for reproducible scientific data analysis. In general these workflows
are designed to be run using Cromwell, although in principal this repo
could also be used to keep track of workflows for any other language or 
system.


### Workflow Folder Organization

Each workflow is stored in its own folder, which contains:
  * The workflow file itself (e.g., *.wdl)
  * An example inputs file (*.json)
  * A readme in Markdown (*.md)

### Data Organization Principals

We subscribe to the philosophy of data organization that each project or
experiment consists of a collection of input files, and that the best way
to analyze that set of files is to create a *manifest* file listing each
of those files. The *workflow* takes that *manifest* as an input, and is
constructed in such a way that it analyzes all of the files within it.

### Validated Tools

A validated tool is one that can be successfully executed using the WDL
and inputs JSON, and is determined by the developer to be ready for use.

#### List of validated tools:

  * Coming soon...

