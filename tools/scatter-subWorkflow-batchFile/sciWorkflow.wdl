workflow sciWorkflow {
    Array[String] chromosomes
    String sampleName

    scatter (chr in chromosomes){
      call printStrings {
        input: 
        sampleName=sampleName,
        chr=chr
      }
    }
    call concatStrings {
        input: 
        files = printStrings.stdout
    }
    output {
      File subworkflow_out = concatStrings.mergedFiles
    }
  }

task printStrings {
    String sampleName
    String chr

    command {
    echo "Analyzing Sample: ${sampleName}, for chromsosome: ${chr}."
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