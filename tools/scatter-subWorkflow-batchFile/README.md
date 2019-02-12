# Scatter over a batchFile but use a subWorkflow that has a scatter as well


The idea here is that one workflow submission would be intended to run one workflow on a arbitrary list of unique samples provided by a batchFile indicating the locations of the input data for each unique sample.  Then the workflow would scatter over each of those samples.  

This works as long as no task in the workflow would benefit from a scatter as well.  In this case any task requiring a scatter would be defined as a subWorkflow and inside that subWorkflow a scatter over an array such as by chromosome or other feature would exist.  It might return an Array[File] or File, depending on the subsequent task in the workflow.  
