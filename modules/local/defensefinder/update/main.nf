process DEFENSEFINDER_UPDATE {
    tag 'single'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/defense-finder:2.0.0--pyhdfd78af_0':
        'biocontainers/defense-finder:2.0.0--pyhdfd78af_0' }"

    output:
    path("defense_db") , emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    defense-finder \\
        update \\
        $args \\
        --models-dir defense_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir defense_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """
}
