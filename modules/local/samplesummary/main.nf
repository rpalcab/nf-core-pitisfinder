process SAMPLESUMMARY {
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
    tuple val(meta), path("${meta.id}_summary.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}_summary.png"), emit: png
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tsv_str = tsv_list.join(',')
    """
    merge_egm.py -s $prefix -g $gbk -t $tsv_str -o ${prefix}_summary
    circos_plot.py -i ${prefix}_summary.gbk -o ${prefix}_summary.png
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
