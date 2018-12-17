import "sciWorkflow.wdl" as sub

workflow batchSubmit {
    File batchFile
    Array[Object] batchInfo = read_objects(batchFile)

  scatter (job in batchInfo){
    call sub.sciWorkflow {
        input: 
        sampleName=job.sampleName, 
        bamLocation=job.bamLocation, 
        bedLocation=job.bedLocation
    }
  }
  output {
      Array[File] batch_output = sciWorkflow.workflow_out
  }
}
