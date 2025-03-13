process INTEGRON_PARSER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
    input:
    tuple val(meta), path(fasta), path(ann), path(res), path(integrons)

    output:
    //tuple val(meta), path("${meta.id}/IS_chr_filtered.tsv"), emit: report
    stdout

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    echo "integron_parser $prefix $fasta $integrons $ann $res"
    integron_parser.py -i $integrons -f $fasta -a $ann -r $res -o .
    """

    stub:
    """
    echo "integron_parser $prefix $integrons $ann $res"
    """
}
