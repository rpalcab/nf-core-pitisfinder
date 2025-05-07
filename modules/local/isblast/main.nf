process IS_BLAST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4':
        'biocontainers/blast:2.16.0--h66d330f_4' }"

    input:
    tuple val(meta), path(fasta)
    path (isdb)

    output:
    tuple val(meta), path("${meta.id}/IS_chr.tsv"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    mkdir ${prefix}
    sed -i 's/>contig/>Chr/g' $fasta
    makeblastdb -in $fasta -dbtype nucl
    blastn -db $fasta -query $isdb \\
        -outfmt "6 qseqid sseqid qstart qend qlen sstart send slen pident qcovhsp length mismatch score evalue" \\
        -evalue 5E-10 -perc_identity 90 -qcov_hsp_perc 90 \\
        -out ${prefix}/IS_chr.tsv

    printf "qseqid\\tsseqid\\tqstart\\tqend\\tqlen\\tsstart\\tsend\\tslen\\tpident\\tqcovhsp\\tlength\\tmismatch\\tscore\\tevalue\\n" > tmp
    cat ${prefix}/IS_chr.tsv >> tmp
    mv tmp ${prefix}/IS_chr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}
    touch ${prefix}/IS_chr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """
}
