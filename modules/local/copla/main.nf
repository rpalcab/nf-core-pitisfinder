process COPLA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:1.0' :
        'docker.io/rpalcab/copla:1.0' }"

    input:
    tuple val(meta), val(plasmid_name), path(file)

    output:
    path "copla/$plasmid_name/", emit: results
    path "copla/$plasmid_name/copla.txt", emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Descargar bases de datos si no existen
    # if [ ! -d "databases" ]; then
    #     bin/download_Copla_databases.sh
    # fi

    mkdir -p copla/${plasmid_name}

    copla \\
         "$file" \\
         /data/app/databases/Copla_RS84/RS84f_sHSBM.pickle \\
         /data/app/databases/Copla_RS84/CoplaDB.fofn \\
         copla/${plasmid_name} > copla/${plasmid_name}/copla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$(echo \$(python3 bin/copla.py --version 2>&1) | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p copla/
    touch copla/copla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$(echo \$(python3 bin/copla.py --version 2>&1) | cut -f2 -d' ')
    END_VERSIONS
    """
}
