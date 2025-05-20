process INTEGRON_PARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
    input:
    tuple val(meta), path(gbk), path(integrons)

    output:
    tuple val(meta), path("${meta.id}/integrons_summary.tsv"), emit: report, optional: true
    tuple val(meta), path("${meta.id}/int_*.fasta"), emit: fasta, optional: true
    tuple val(meta), path("${meta.id}/int_*.tsv"), emit: tsv, optional: true
    tuple val(meta), path("${meta.id}/int_*.gbk"), emit: gbk, optional: true
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    integron_parser.py -i $integrons -a $gbk -o ${prefix}
    """

    stub:
    """
    touch integrons_summary.tsv
    """
}
