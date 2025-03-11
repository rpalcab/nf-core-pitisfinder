process IS_PARSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://continuumio/anaconda3:latest':
        'docker.io/continuumio/anaconda3:latest' }"
    input:
    tuple val(meta), path(report_raw)

    output:
    tuple val(meta), path("${meta.id}/IS_chr_filtered.tsv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    mkdir ${prefix}
    IS_parser.py -i $report_raw -o ${prefix}
    """

    stub:
    """
    mkdir ${prefix}
    touch ${prefix}/IS_chr_filtered.tsv
    """
}
