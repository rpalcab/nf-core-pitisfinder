process VISUALIZE_LINEAR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), val(mge), path(gbk)
    val(outdir)

    output:
    tuple val(meta), val(mge), path("${outdir}/${mge.id}.png"), emit: png
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $outdir
    linear_plot.py -i $gbk -o $outdir/${mge.id}.png
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${mge.id}.png
    """
}
