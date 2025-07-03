process PLASMIDMARKERS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(biomarkers_mobsuite), path(biomarkers_copla), path(gbk)

    output:
    tuple val(meta), path("${gbk.baseName}_pl.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    def outfile = "${gbk.baseName}_pl.gbk"
    """
    echo -e "Contig\\tstart\\tend\\tstrand\\tGene\\tIdentity\\tCoverage\\ti-evalue\\ttag" > merged_biomarkers_copla.tsv
    tail -q -n+2 $biomarkers_copla >> merged_biomarkers_copla.tsv
    merge_plasmid_ann.py -g $gbk -m $biomarkers_mobsuite -c merged_biomarkers_copla.tsv -o "$outfile"
    """
}
