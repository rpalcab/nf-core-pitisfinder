process SAMPLESUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), path(tsv_list), path(gbk_list), path(gbk)

    output:
    tuple val(meta), path("${meta.id}_markers.gbk"), emit: gbk, optional: true
    tuple val(meta), path("${meta.id}_summary.tsv"), emit: tsv, optional: true
    tuple val(meta), path("${meta.id}_summary.png"), emit: png, optional: true
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tsv_str = tsv_list.join(',')
    def tsv_opt = tsv_str ? "-t ${tsv_str}" : ""
    def gbk_merged = tsv_str ? "${prefix}_summary.gbk" : gbk
    def gbk_str = gbk_list.join(',')
    def gbk_opt = "-u ${gbk_str}"
    """
    merge_egm.py \\
            -s $prefix \\
            -g $gbk \\
            ${tsv_opt} \\
            -o ${prefix}_summary

    merge_gbks.py -g ${prefix}_summary.gbk $gbk_opt -o ${prefix}_markers.gbk

    circos_plot.py -i ${prefix}_markers.gbk -o ${prefix}_summary.png
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_summary.png
    touch ${prefix}_summary.tsv
    touch ${prefix}_summary.gbk
    """
}
