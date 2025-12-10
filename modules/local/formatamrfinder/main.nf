process FORMATAMRFINDER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${meta.id}_reformated.tsv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    format_amrfinder.py \\
        $args \\
        -i $tsv \\
        -o ${prefix}_reformated.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_reformated.tsv
    """
}
