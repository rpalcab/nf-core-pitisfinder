process PLASMIDMARKERS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(biomarkers), path(gbk)

    output:
    tuple val(meta), path("${gbk.baseName}_pl.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    def outfile = "${gbk.baseName}_pl.gbk"
    """
    merge_plasmid_ann.py -g $gbk -b $biomarkers -o "$outfile"
    """
}
