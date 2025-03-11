process IS_BLAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/ncbi/blast:2.16.0':
        'docker.io/ncbi/blast:2.16.0' }"
    input:
    tuple val(meta), path(fasta)
    path (isdb)

    output:
    tuple val(meta), path("${meta.id}/IS_chr_raw.tsv"), emit: report
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
        -out ${prefix}/IS_chr_raw.tsv
    sed -i '1s/^/qseqid\\tsseqid\\tqstart\\tqend\\tqlen\\tsstart\\tsend\\tslen\\tpident\\tqcovhsp\\tlength\\tmismatch\\tscore\\tevalue\\n/' ${prefix}/IS_chr_raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}
    touch ${prefix}/IS_chr_raw.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """
}
