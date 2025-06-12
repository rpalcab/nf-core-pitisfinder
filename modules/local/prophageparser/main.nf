process PROPHAGEPARSER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"

    input:
    tuple val(meta), path(genomic_gbk), path(provirus), path(taxonomy)

    output:
    tuple val(meta), path("prophage_summary.tsv"), emit: summary
    tuple val(meta), path("phage_*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("phage_*.tsv"), emit: tsv, optional: true
    tuple val(meta), path("phage_*.gbk"), emit: gbk, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prophage_parser.py -p $provirus -t $taxonomy -a $genomic_gbk -o .
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch prophage_summary.tsv
    """
}
