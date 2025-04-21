process ICEBERG_ICESEARCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4':
        'biocontainers/blast:2.16.0--h66d330f_4' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}/blast.tsv"), emit: tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeblastdb -in $fasta -dbtype nucl

    mkdir ${prefix}
    blastn -db $fasta -query $db \\
        -outfmt "6 qseqid sseqid qstart qend qlen sstart send slen pident qcovhsp length mismatch score evalue" \\
        -evalue 5E-10 -perc_identity 80 -qcov_hsp_perc 60 \\
        -out ${prefix}/blast.tsv
    sed -i '1s/^/qseqid\\tsseqid\\tqstart\\tqend\\tqlen\\tsstart\\tsend\\tslen\\tpident\\tqcovhsp\\tlength\\tmismatch\\tscore\\tevalue\\n/' ${prefix}/blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iceberg: \$(blastn -version | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """
}
