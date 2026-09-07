process PROPHAGEMARKERS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"

    input:
    tuple val(meta), path(biomarkers), path(provirus), path(gbk)

    output:
    tuple val(meta), path("${gbk.baseName}_ph.gbk"), emit: gbk, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile = "${gbk.baseName}_ph.gbk"
    """
    merge_phage_ann.py -g $gbk -b $biomarkers -i $provirus -o "$outfile"
    """
}
