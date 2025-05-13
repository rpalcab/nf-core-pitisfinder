process PHIGARO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phigaro:2.4.0--pyhdfd78af_0':
        'biocontainers/phigaro:2.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)
    path(config)

    output:
    path("${meta.id}/${meta.id}.phigaro.bed")  , optional:true, emit: bed
    path("${meta.id}/${meta.id}.phigaro.fasta"), optional:true, emit: fasta
    path("${meta.id}/${meta.id}.phigaro.gff3") , optional:true, emit: gff3
    path("${meta.id}/${meta.id}.phigaro.html") , optional:true, emit: html
    path("${meta.id}/${meta.id}.phigaro.tsv")  , optional:true, emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    phigaro \\
        $args \\
        -f $fasta \\
        -o $prefix \\
        -c $config \\
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phigaro: \$(phigaro -V  |& sed '1!d ; s/phigaro //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phigaro: \$(phigaro -V  |& sed '1!d ; s/phigaro //')
    END_VERSIONS
    """
}