process VISUALIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), val(plasmid_name), path(gbk)

    output:
    tuple val(meta), val(plasmid_name), path("*.png"), emit: bam
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circos_plot.py -i $gbk -o ${plasmid_name}.png
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${plasmid_name}.png
    """
}
