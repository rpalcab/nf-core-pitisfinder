process GENOMAD_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.11.0--pyhdfd78af_0':
        'biocontainers/genomad:1.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    path  genomad_db

    output:
    tuple val(meta), path("*_annotate/*_genes.tsv")       , emit: gene_annotation
    tuple val(meta), path("*_annotate/*_proteins.faa")    , emit: faa
    tuple val(meta), path("*_annotate")                   , emit: outdir
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        annotate \\
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
    mkdir ${prefix}_annotate
    touch ${prefix}_annotate/*_genes.tsv
    touch ${prefix}_annotate/*_proteins.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
