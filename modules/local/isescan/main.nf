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
    tuple val(meta), path("*.tsv")   , emit: tsv, optional: true
    tuple val(meta), path("*.gff")   , emit: gff, optional: true
    tuple val(meta), path("*.is.fna"), emit: isfna, optional: true
    tuple val(meta), path("*.orf.fna"), emit: orffna, optional: true
    tuple val(meta), path("*.orf.faa"), emit: orffaa, optional: true
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
        --output . \\
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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$(isescan.py --version |& sed '1!d ; s/isescan //')
    END_VERSIONS
    """
}
