workflow sciWorkflow {
    String sampleName 
    File bamLocation
    File bedLocation
    call testTask {
        input: 
        in1=sampleName, 
        in2=bamLocation, 
        in3=bedLocation
    }
    output {
        File workflow_out = testTask.testTask_out
    }
  }

task testTask {
    String in1
    String in2
    String in3
    command {
    echo ${in1}
    echo ${in2}
    echo ${in3}
    }
    output {
        File testTask_out = stdout()
    }
}