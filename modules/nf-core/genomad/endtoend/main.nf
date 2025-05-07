process GENOMAD_ENDTOEND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.11.0--pyhdfd78af_0':
        'biocontainers/genomad:1.11.0--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    path  genomad_db

    output:
    tuple val(meta), path("${meta.id}/*_aggregated_classification/*_aggregated_classification.tsv")    , emit: aggregated_classification   , optional: true
    tuple val(meta), path("${meta.id}/*_annotate/*_taxonomy.tsv")                                      , emit: taxonomy
    tuple val(meta), path("${meta.id}/*_find_proviruses/*_provirus.tsv")                               , emit: provirus
    tuple val(meta), path("${meta.id}/*_score_calibration/*_compositions.tsv")                         , emit: compositions                , optional: true
    tuple val(meta), path("${meta.id}/*_score_calibration/*_calibrated_aggregated_classification.tsv") , emit: calibrated_classification   , optional: true
    tuple val(meta), path("${meta.id}/*_summary/*_plasmid.fna.gz")                                     , emit: plasmid_fasta
    tuple val(meta), path("${meta.id}/*_summary/*_plasmid_genes.tsv")                                  , emit: plasmid_genes
    tuple val(meta), path("${meta.id}/*_summary/*_plasmid_proteins.faa.gz")                            , emit: plasmid_proteins
    tuple val(meta), path("${meta.id}/*_summary/*_plasmid_summary.tsv")                                , emit: plasmid_summary
    tuple val(meta), path("${meta.id}/*_summary/*_virus.fna.gz")                                       , emit: virus_fasta
    tuple val(meta), path("${meta.id}/*_summary/*_virus_genes.tsv")                                    , emit: virus_genes
    tuple val(meta), path("${meta.id}/*_summary/*_virus_proteins.faa.gz")                              , emit: virus_proteins
    tuple val(meta), path("${meta.id}/*_summary/*_virus_summary.tsv")                                  , emit: virus_summary
    path "versions.yml"                                                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        end-to-end \\
        $fasta \\
        $prefix \\
        $genomad_db \\
        --threads $task.cpus \\
        $args

    gzip $prefix/**/*.fna
    gzip $prefix/**/*.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    """
    mkdir $prefix/${filename}_aggregated_classification
    touch $prefix/${filename}_aggregated_classification/${filename}_aggregated_classification.tsv
    mkdir $prefix/${filename}_annotate
    touch $prefix/${filename}_annotate/${filename}_taxonomy.tsv
    mkdir $prefix/${filename}_find_proviruses
    touch $prefix/${filename}_find_proviruses/${filename}_provirus.tsv
    mkdir $prefix/${filename}_marker_classification
    mkdir $prefix/${filename}_nn_classification
    mkdir $prefix/${filename}_score_calibration
    touch $prefix/${filename}_score_calibration/${filename}_calibrated_aggregated_classification.tsv
    touch $prefix/${filename}_score_calibration/${filename}_compositions.tsv
    mkdir $prefix/${filename}_summary
    touch $prefix/${filename}_summary/${filename}_plasmid.fna.gz
    touch $prefix/${filename}_summary/${filename}_plasmid_genes.tsv
    touch $prefix/${filename}_summary/${filename}_plasmid_proteins.faa.gz
    touch $prefix/${filename}_summary/${filename}_plasmid_summary.tsv
    touch $prefix/${filename}_summary/${filename}_virus.fna.gz
    touch $prefix/${filename}_summary/${filename}_virus_genes.tsv
    touch $prefix/${filename}_summary/${filename}_virus_proteins.faa.gz
    touch $prefix/${filename}_summary/${filename}_virus_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
