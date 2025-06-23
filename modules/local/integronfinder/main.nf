process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/integron_finder:2.0.5--pyhdfd78af_0':
        'biocontainers/integron_finder:2.0.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.integrons"), emit: integrons
    tuple val(meta), path("${meta.id}.summary")  , emit: summary
    tuple val(meta), path("integron_finder.out") , emit: log         , optional: true
    tuple val(meta), path("*.gbk")               , emit: gbk         , optional: true
    tuple val(meta), path("*.pdf")               , emit: pdf         , optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    integron_finder \\
                $fasta \\
                $args \\
                --cpu $task.cpus \\
                --outdir .
    mv Results_Integron_Finder_${prefix}/* .
    rmdir Results_Integron_Finder_${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(integron_finder --version | head -n1 | sed 's/^integron_finder version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.integrons
    touch ${prefix}.summary"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(integron_finder --version | head -n1 | sed 's/^integron_finder version //')
    END_VERSIONS
    """
}
