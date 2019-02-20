
## The goal here would be to scatter over a batch file in S3 of 400+ independent datasets,
## then perform using the same workflow options and paramters, a series of tasks, which could
## include a scatter/scatter-gather format.  Independent jobs should complete if possible 
## regardless of other dataset-scatter shards.  
workflow sciWorkflow {
    File batchFile
    Array[String] chromosomes
    Array[Object] batchInfo = read_objects(batchFile)

  scatter (job in batchInfo){
    String sampleName = job.sampleName
    String metadata1 = job.metadata1

    scatter (chr in chromosomes){
      call printStrings {
        input:
        sample = sampleName,
        metadata = metadata1,
        chr = chr
      }
    }
    call concatStrings{
      input: 
      files = printStrings.stdout
    }
  }
  output {
    Array[File] doublescatter = concatStrings.mergedFiles
  }
}

task printStrings {
    String sample
    String metadata
    String chr

  command {
    echo "Analyzing Sample: ${sample}, with metadata ${metadata}, for chromosome: ${chr}."
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

task concatStrings {
    Array[File] files

  command {
    cat ${write_lines(files)} > merged.txt
    }
  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: "1"
  }
  output {
    File mergedFiles = "merged.txt"
  }

}