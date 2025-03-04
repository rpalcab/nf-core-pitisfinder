process COPLA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:latest' :
        'docker.io/rpalcab/copla:latest' }"

    input:
    tuple val(meta), path(plasmid_file)

    output:
    path "copla/${meta.id}/", emit: results
    path "copla.txt", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Descargar bases de datos si no existen
    # if [ ! -d "databases" ]; then
    #     bin/download_Copla_databases.sh
    # fi

    python3 bin/copla.py \\
         "$plasmid_file" \\
         databases/Copla_RS84/RS84f_sHSBM.pickle \\
         databases/Copla_RS84/CoplaDB.fofn \\
         copla/${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$(echo \$(python3 bin/copla.py --version 2>&1) | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$(echo \$(python3 bin/copla.py --version 2>&1) | cut -f2 -d' ')
    END_VERSIONS
    """
}
