workflow myWorkflow {
    call myTask
}

task myTask {
    String name
    command {
        echo "hello ${name}"
    }
    output {
        File response = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
