process COPLA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
<<<<<<< Updated upstream
        'docker://rpalcab/copla:1.0' :
        'docker.io/rpalcab/copla:1.0' }"
=======
        'docker://rpalcab/copla:1.0light' :
        'docker.io/rpalcab/copla:1.0light' }"
>>>>>>> Stashed changes

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
    download_Copla_databases
    # post_install_test

    mkdir -p copla/${plasmid_name}

    copla \\
         "$file" \\
         databases/Copla_RS84/RS84f_sHSBM.pickle \\
         databases/Copla_RS84/CoplaDB.fofn \\
         copla/${plasmid_name} > copla/${plasmid_name}/copla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$( copla --version | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p copla/
    touch copla/copla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$( copla --version | cut -f2 -d' ')
    END_VERSIONS
    """
}
