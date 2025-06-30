process MERGE_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

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
