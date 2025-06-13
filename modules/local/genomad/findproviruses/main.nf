process GENOMAD_FINDPROVIRUSES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.11.0--pyhdfd78af_0':
        'biocontainers/genomad:1.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta), path(annotate_dir)
    path  genomad_db

    output:
    tuple val(meta), path("*_find_proviruses")                        , emit: outdir
    tuple val(meta), path("*_find_proviruses/*_provirus.tsv")         , emit: provirus
    tuple val(meta), path("*_find_proviruses/*_provirus_genes.tsv")   , emit: provirus_genes
    tuple val(meta), path("*_find_proviruses/*_provirus_taxonomy.tsv"), emit: taxonomy
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        find-proviruses \\
        $fasta \\
        . \\
        $genomad_db \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_find_proviruses
    touch ${prefix}_find_proviruses/${prefix}_provirus.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
