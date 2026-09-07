process INTEGRONMARKERS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"

    input:
    tuple val(meta), path(biomarkers), path(gbk)

    output:
    tuple val(meta), path("${gbk.baseName}_int.gbk"), emit: gbk, optional: true

    script:
    def prefix = "${meta.id}"
    def outfile = "${gbk.baseName}_int.gbk"
    """
    merge_integron_ann.py -g $gbk -b $biomarkers -o "$outfile"
    """
}
