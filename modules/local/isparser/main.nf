process IS_PARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
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
