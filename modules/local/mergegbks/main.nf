process MERGEGBKS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), path(summary_gbk), path(gbk_list)

    output:
    tuple val(meta), path("${meta.id}_markers.gbk"), emit: gbk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gbk_str = gbk_list.join(',')
    def gbk_opt = gbk_str ? "-u ${gbk_str}" : ""
    """
    merge_gbks.py -g ${summary_gbk} $gbk_opt -o ${prefix}_markers.gbk
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_markers.gbk
    """
}
