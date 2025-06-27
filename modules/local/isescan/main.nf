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
    tuple val(meta), path("${fasta}.tsv")        , emit: tsv, optional: true
    tuple val(meta), path("${fasta}.gff")        , emit: gff, optional: true
    tuple val(meta), path("${fasta}.is.fna")     , emit: isfna, optional: true
    tuple val(meta), path("${fasta}.orf.fna")    , emit: orffna, optional: true
    tuple val(meta), path("${fasta}.orf.faa")    , emit: orffaa, optional: true
    tuple val(meta), path("IS_summary.tsv")      , emit: summary, optional: true
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

    if [ -f ${fasta}.tsv ]; then
            cut -f 1,3-6 ${fasta}.tsv | awk -F'\\t' 'BEGIN {OFS="\\t"} NR==1 {print \$0, "AMR\\tVF"; next} {print \$0, "\\t"}' | sed -e 's/seqID/Contig/' -e 's/cluster/Name/' -e 's/isBegin/Start/' -e 's/isEnd/End/' -e 's/isLen/Length/'> IS_summary.tsv
    fi

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
