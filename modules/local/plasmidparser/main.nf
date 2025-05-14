process PLASMID_PARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"
    input:
    tuple val(meta), val(plasmid_name), path(qry_info), path(ptu), path(gff), path(amr), path(mobsuite)

    output:
    tuple val(meta), path("$plasmid_name/plasmids_summary.tsv"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    plasmid_parser.py -m $mobsuite -q $qry_info -p $ptu -a $gff -r $amr -n $plasmid_name -o $plasmid_name
    """

    stub:
    """
    mkdir ${prefix}
    touch ${prefix}/plasmids_summary.tsv
    """
}
