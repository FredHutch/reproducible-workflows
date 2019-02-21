
## The goal here would be to scatter over a batch file in S3 of 400+ independent datasets,
## then perform using the same workflow options and paramters, a series of tasks, which could
## include a scatter/scatter-gather format.  Independent jobs should complete if possible 
## regardless of other dataset-scatter shards.  
workflow batchTesting {
    File batchFile
    Array[String] chromosomes
    Array[Object] batchInfo = read_objects(batchFile)

  scatter (job in batchInfo){
    String sampleName = job.sampleName
    String metadata1 = job.metadata1

      call printStrings {
        input:
        sample = sampleName,
        metadata = metadata1,
        chr = chromosomes
      }

  }
  output {
    Array[File] singlescatter = printStrings.stdout
  }
}

task printStrings {
    String sample
    String metadata
    Array[String] chr

  command {
    echo "Analyzing Sample: ${sample}, with metadata ${metadata}, for chromosome: ${sep=" " chr}."
    }
  
  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: "1"
  }
  output {
    File stdout=stdout()
  }

}
