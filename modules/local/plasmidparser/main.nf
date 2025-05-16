process PLASMID_PARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
    input:
    tuple val(meta), val(plasmid_name), path(qry_info), path(ptu), path(gbk), path(mobsuite), path(contig_report)

    output:
    tuple val(meta), path("${plasmid_name}.tsv"), emit: report
    tuple val(meta), path("${plasmid_name}.gbk"), emit: gbk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    plasmid_parser.py -m $mobsuite -r $contig_report -q $qry_info -p $ptu -g $gbk -n $plasmid_name -o .
    """

    stub:
    """
    mkdir ${prefix}
    touch ${prefix}/${plasmid_name}.tsv
    touch ${prefix}/${plasmid_name}.gbk
    """
}
