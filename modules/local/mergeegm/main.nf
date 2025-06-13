process MERGEEGM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), path(tsv_list), path(gbk)

    output:
    tuple val(meta), path("${meta.id}_summary.gbk"), emit: gbk
    tuple val(meta), path("${meta.id}_summary.tsv"), emit: tsv, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tsv_str = tsv_list.join(',')
    def tsv_opt = tsv_str ? "-t ${tsv_str}" : ""
    """
    merge_egm.py \\
        -s $prefix \\
        -g $gbk \\
        ${tsv_opt} \\
        -o ${prefix}_summary
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.gbk
    touch ${prefix}_summary.tsv
    """
}
