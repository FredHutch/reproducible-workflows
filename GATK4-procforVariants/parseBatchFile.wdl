workflow parseBatchFile {
  File batchfile


call findJobs {
    input:
      docker_image = gotc_docker,
      bwa_path = gotc_path,
      preemptible_tries = preemptible_tries
  }


task findJobs {
  command <<<
    python <<CODE
    print('\t'.join(["key_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    CODE
  >>>
  output {
    Array[Object] my_obj = read_objects(stdout())
  }
}


workflow main {
  call findJobs {}
  scatter (individualJob in jobList){
  Array[String] theseParams = individualJob 
  call breakout {}
  }



task findJobs {
  File batchfile
  command <<<
    python <<CODE
    print('\t'.join(["key_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    print('\t'.join(["value_{}".format(i) for i in range(3)]))
    CODE
  >>>
  output {
    Array[Object] jobList = read_objects(stdout())
  }
}

task breakout {
  command <<<
    echo ${write_objects(jobList)}
  >>>
  output {
    Array[String] jobParams = read_objects(stdout())}
  }

}
