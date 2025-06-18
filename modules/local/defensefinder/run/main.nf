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
    cat > test.py<<EOF
    import os
    file_lastver = os.path.join(os.environ["HOME"], ".defensefinder_model_lastversion")
    print(file_lastver)
    print(os.environ["HOME"])
    EOF

    python3 test.py
    find . -name ".defensefinder_model_lastversion"

    defense-finder \\
        run \\
        $args \\
        -w $task.cpus \\
        --models-dir $db \\
        $fasta

    grep '^>' ${prefix}.prt | awk -F ' # ' '{gsub(/^>/, "", \$1); print \$1, \$2, \$3, \$4}' > gene_positions.tsv

    awk '
        BEGIN {
            FS = OFS = "\\t"
        }
        FNR==NR {
            print \$0
            pos_start[\$1] = \$2
            pos_end[\$1] = \$3
            strand[\$1] = \$4
            next
        }
        FNR==1 {
            print \$0, "start_pos", "end_pos", "strand"
            next
        }
        {
            id = \$2
            print \$0, (id in pos_start ? pos_start[id] : "NA"), (id in pos_end ? pos_end[id] : "NA"), (id in strand ? strand[id] : "NA")
        }
    ' gene_positions.tsv ${prefix}_defense_finder_genes.tsv > ${prefix}_defense_finder_genes_pos.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_defense_finder_genes.tsv
    touch ${prefix}_defense_finder_hmmer.tsv
    touch ${prefix}_defense_finder_systems.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defensefinder: \$(defense-finder version |& sed 's/Using DefenseFinder version //')
    END_VERSIONS
    """
}
