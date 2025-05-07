process COPLA_COPLADBDOWNLOAD {
    tag "copladb"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:1.0light' :
        'docker.io/rpalcab/copla:1.0light' }"

    output:
    path "databases", emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Downloading Copla DB"
    download_Copla_databases
    """

    stub:
    """
    mkdir -p databases/
    """
}
