process MACSYFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macsyfinder:2.1.4--pyhdfd78af_1':
        'biocontainers/macsyfinder:2.1.4--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(faa)
    tuple val(model), path(model_path)

    output:
    tuple val(meta), path("all_best_solutions.tsv")  , emit: best
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    macsyfinder \\
        $args \\
        --sequence-db $faa \\
        -m $model all \\
        --models-dir $model_path \\
        -o tmp/ \\
        -w $task.cpus
    mv tmp/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch all_best_solutions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """
}
