process DEFENSEFINDER_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/defense-finder:2.0.1--pyhdfd78af_0':
        'biocontainers/defense-finder:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_defense_finder_genes_pos.tsv")  , emit: genes
    tuple val(meta), path("${meta.id}_defense_finder_hmmer.tsv")  , emit: hmmer
    tuple val(meta), path("${meta.id}_defense_finder_systems.tsv"), emit: systems
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    defense-finder \\
        run \\
        $args \\
        -w $task.cpus \\
        --models-dir $db \\
        $fasta

    tail -n +2 ${prefix}_defense_finder_genes.tsv | sort -k2,2 > ${prefix}_defense_finder_genes_sorted.tsv
    grep '^>' ${prefix}.prt  | awk -F ' # ' '{gsub(/^>/, "", \$1); print \$1, \$2, \$3, \$4}' | sed 's/\\s/\\t/g' | sort -k1,1 > gene_positions_sorted.tsv
    echo -e "replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tlocus_num\\tsys_wholeness\\tsys_score\\tsys_occ\\thit_gene_ref\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov\\thit_seq_cov\\thit_begin_match\\thit_end_match\\tcounterpart\\tused_in\\ttype\\tsubtype\\tactivity\\tstart\\tend\\tframe" > ${prefix}_defense_finder_genes_pos.tsv
    awk -F'\t' 'NR==FNR {pos[\$1]=\$0; next}
                \$2 in pos {
                    split(pos[\$2], a, "\\t");
                    print \$0 "\\t" a[2] "\\t" a[3] "\\t" a[4]
                }' gene_positions_sorted.tsv ${prefix}_defense_finder_genes_sorted.tsv >> ${prefix}_defense_finder_genes_pos.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_defense_finder_genes_pos.tsv
    touch ${prefix}_defense_finder_hmmer.tsv
    touch ${prefix}_defense_finder_systems.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """
}
