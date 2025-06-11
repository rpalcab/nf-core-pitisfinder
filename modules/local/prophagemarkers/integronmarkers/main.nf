process INTEGRONMARKERS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(biomarkers), path(gbk)

    output:
    tuple val(meta), path("${gbk.baseName}_int.gbk"), emit: gbk, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def outfile = "${gbk.baseName}_int.gbk"
    """
    merge_integron_ann.py -g $gbk -b $biomarkers -o "$outfile"
    """
}
