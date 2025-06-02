process MERGE_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(amr), path(gbk)

    output:
    tuple val(meta), path("${meta.id}_merged.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    """
    merge_amr.py -t $amr -g $gbk -o "$prefix"_merged.gbk
    """
}
