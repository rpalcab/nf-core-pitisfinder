process FILTER_FAA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(faa), path(gbk), path(chr)

    output:
    tuple val(meta), path("${meta.id}/*_filtered.faa"), emit: chr

    script:
    def prefix = "${meta.id}"
    """
    cut -f2,5 $chr | grep "chromosome" | cut -f2 > chr.txt
    faa_contig_filter.py -f $faa -g $gbk -c chr.txt -o $prefix
    """
}