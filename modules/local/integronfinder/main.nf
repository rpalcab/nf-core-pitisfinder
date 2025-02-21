process INTEGRON_FINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/integron_finder:latest':
        'docker.io/rpalcab/integron_finder:latest' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}"), emit: outdir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    echo 'integron_finder'
    integron_finder $fasta --cpu ${task.cpus} --outdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(echo \$(integron_finder --version 2>&1) | head -n1 | sed 's/^integron_finder version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(echo \$(integron_finder --version 2>&1) | head -n1 | sed 's/^integron_finder version //')
    END_VERSIONS
    """
}
