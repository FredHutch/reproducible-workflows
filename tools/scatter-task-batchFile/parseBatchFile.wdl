workflow parseBatchFile {
  File batchFile
  Array[Object] batchInfo = read_objects(batchFile)
  scatter (job in batchInfo){
    String sampleName = job.sampleName
    File bamLocation = job.bamLocation
    File bedLocation = job.bedLocation
    call test {
        input: in1=sampleName, in2=bamLocation, in3=bedLocation
    }
  }
}
# echo ${sep=' ' in}
task test {
    String in1
    String in2
    String in3
    command {
    echo ${in1}
    echo ${in2}
    echo ${in3}
    }
    output {
        File item_out = stdout()
    }
}