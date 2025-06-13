process DEFENSEFINDER_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/defense-finder:2.0.0--pyhdfd78af_0':
        'biocontainers/defense-finder:2.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_defense_finder_genes.tsv"), emit: genes
    tuple val(meta), path("${meta.id}_defense_finder_hmmer.tsv"), emit: hmmer
    tuple val(meta), path("${meta.id}_defense_finder_systems.tsv"), emit: systems
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > test.py<<EOF
    import os
    file_lastver = os.path.join(os.environ["HOME"], ".defensefinder_model_lastversion")
    print(file_lastver)
    print(os.environ["HOME"])
    EOF

    python3 test.py
    find . -name ".defensefinder_model_lastversion"

    defense-finder \\
        run \\
        $args \\
        -w $task.cpus \\
        --models-dir $db \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_defense_finder_genes.tsv
    touch ${prefix}_defense_finder_hmmer.tsv
    touch ${prefix}_defense_finder_systems.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """
}
