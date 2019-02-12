import "sciWorkflow.wdl" as sub1

workflow batchSubmit {
    Array[String] batchString

  scatter (job in batchString){
    Array[String] chromosomes = ["chr1", "chr2"]

    call sub1.sciWorkflow {
        input: 
        sampleName=job,
        chromosomes=chromosomes
    }
  }
  output {
      Array[File] batch_output = sciWorkflow.subworkflow_out
  }
}
