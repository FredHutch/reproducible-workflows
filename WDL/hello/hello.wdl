workflow myWorkflow {
    call myTask
}

task myTask {
    string name
    command {
        echo "hello ${name}"
    }
    output {
        File response = stdout()
    }
}
