process PHASTEST_PHASTESTDBDOWNLOAD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/phastest_single:1.1' :
        'docker.io/rpalcab/phastest_single:1.1' }"

    output:
    path "DB/", emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget https://phastest.ca/download_file/docker-database.zip
    unzip docker-database.zip
    """

    stub:
    """
    mkdir -p DB/
    """
}
