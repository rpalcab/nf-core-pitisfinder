process ISESCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isescan:1.7.2.3--h7b50bb2_3':
        'biocontainers/isescan:1.7.2.3--h7b50bb2_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.tsv")   , emit: tsv
    tuple val(meta), path("${meta.id}/*.gff")   , emit: gff
    tuple val(meta), path("${meta.id}/*.is.fna"), emit: isfna
    tuple val(meta), path("${meta.id}/*.orf.fna"), emit: orffna
    tuple val(meta), path("${meta.id}/*.orf.faa"), emit: orffaa
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isescan.py \\
        $args \\
        --nthread $task.cpus \\
        --output $prefix \\
        --seqfile $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$(isescan.py --version |& sed '1!d ; s/isescan //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$(isescan.py --version |& sed '1!d ; s/isescan //')
    END_VERSIONS
    """
}
