process ICEFINDER2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/iceberg2:1.0' :
        'docker.io/rpalcab/iceberg2:1.0' }"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("${meta.id}"), emit: result
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARNING: Ensure to update software version in this line if you modify the container/environment.
    def icefinder2_version = "1.0"
    """
    ls
    ./ICEfinder2.py \\
        -i $gbk \\
        -t Single

    mv result ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        icefinder2: \$(echo "${icefinder2_version}")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        icefinder2: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
