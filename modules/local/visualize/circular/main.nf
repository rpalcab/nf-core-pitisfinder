process VISUALIZE_CIRCULAR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), val(name), path(gbk), val(mobsuite_report)

    output:
    tuple val(meta), val(name), path("*.png"), emit: png

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = mobsuite_report ? "-r ${mobsuite_report}" : ""
    """
    circos_plot.py -i $gbk $args -o ${name}.png
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${name}.png
    """
}
