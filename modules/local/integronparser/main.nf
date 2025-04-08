process INTEGRON_PARSER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
    input:
    tuple val(meta), path(fasta), path(gene_ann), path(prot_ann), path(res), path(integrons)

    output:
    tuple val(meta), path("integrons_summary.tsv"), emit: report, optional: true
    tuple val(meta), path("int_*.fasta"), emit: fastas, optional: true
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    integron_parser.py -i $integrons -f $fasta -a $gene_ann -r $res -s ${meta.id} -o .
    """

    stub:
    """
    touch integrons_summary.tsv
    """
}
