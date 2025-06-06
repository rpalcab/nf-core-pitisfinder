process PLASMID_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reports), path(gbks_plasmids), path(gbk_genome)

    output:
    tuple val(meta), path("plasmid_summary.tsv"), emit: summary

    script:
    def prefix = "${meta.id}"
    """
    echo -e "Contig\tName\tStart\tEnd\tAMR" > plasmid_summary.tsv
    tail -q -n+2 $reports | cut -f2,6,3,4,19 | awk -F'\t' 'BEGIN{OFS="\t"} {print \$1, \$2":"\$3, 1, \$4, \$5}' >> plasmid_summary.tsv
    """
}
