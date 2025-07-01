process GETCOORDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:1.0light' :
        'docker.io/rpalcab/copla:1.0light' }"

    input:
    tuple val(meta), val(plasmid_name), path(faa), val(conj), path(mob, stageAs: "mob_results.tsv"), path(rep, stageAs: "rep_results.tsv")

    output:
    tuple val(meta), path("${plasmid_name}_copla_biomarkers.tsv"), emit: biomarkers

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def conj_param = conj ? "-c $conj" : ""
    """
    echo -e "gene_id\\tstart\\tend\\tstrand" > gene_positions.tsv
    grep '^>' ${faa} | awk -F ' # ' '{gsub(/^>/, "", \$1); print \$1, \$2, \$3, \$4}' | sed 's/\\s/\\t/g' >> gene_positions.tsv
    copla_coords.py -p gene_positions.tsv -m $mob -r $rep  $conj_param -o ${plasmid_name}_copla_biomarkers.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_copla_biomarkers.tsv
    """
}
