process PLASMIDFINDER {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidfinder:2.1.6--py310hdfd78af_1':
        'biocontainers/plasmidfinder:2.1.6--py310hdfd78af_1' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("${meta.id}/*.json")                 , emit: json
    tuple val(meta), path("${meta.id}/*.txt")                  , emit: txt
    tuple val(meta), path("${meta.id}/*.tsv")                  , emit: tsv
    tuple val(meta), path("${meta.id}/*-hit_in_genome_seq.fsa"), emit: genome_seq
    tuple val(meta), path("${meta.id}/*-plasmid_seqs.fsa")     , emit: plasmid_seq
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.1.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir $prefix
    plasmidfinder.py \\
        $args \\
        -i $seqs \\
        -o $prefix \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv ${prefix}/data.json ${prefix}/${prefix}.json
    mv ${prefix}/results.txt ${prefix}/${prefix}.txt
    mv ${prefix}/results_tab.tsv ${prefix}/${prefix}.tsv
    mv ${prefix}/Hit_in_genome_seq.fsa ${prefix}/${prefix}-hit_in_genome_seq.fsa
    mv ${prefix}/Plasmid_seqs.fsa ${prefix}/${prefix}-plasmid_seqs.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: $VERSION
    END_VERSIONS
    """
}
