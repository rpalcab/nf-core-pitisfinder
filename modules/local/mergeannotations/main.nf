process MERGE_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/pitis_parser:1.0':
        'docker.io/rpalcab/pitis_parser:1.0' }"

    input:
    tuple val(meta), path(gbk), path(amr), path(vf), path(df)

    output:
    tuple val(meta), path("${meta.id}_merged.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    """
    merge_annotation.py -g $gbk -a $amr -v $vf -d $df -o "$prefix"_merged.gbk
    """
}
