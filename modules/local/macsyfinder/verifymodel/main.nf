process VERIFYMODEL {
    tag "verifymodel"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macsyfinder:2.1.4--pyhdfd78af_1':
        'biocontainers/macsyfinder:2.1.4--pyhdfd78af_1' }"

    input:
    val(model)

    output:
    tuple val(model), path("models/")           , emit: model
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def dwnl_model = model.contains("/") ? model.tokenize("/")[0] : model
    """
    macsydata install ${dwnl_model} -t models/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir models/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsyfinder: \$(macsyfinder --version | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """
}
