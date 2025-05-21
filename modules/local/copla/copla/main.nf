process COPLA_COPLA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:1.0light' :
        'docker.io/rpalcab/copla:1.0light' }"

    input:
    tuple val(meta), val(plasmid_name), path(file)
    path (copladb)

    output:
    // tuple val(meta), val(plasmid_name), path("$plasmid_name/"), emit: results
    tuple val(meta), val(plasmid_name), path("$plasmid_name/*.qry_info.tsv"), emit: query
    tuple val(meta), val(plasmid_name), path("$plasmid_name/*.ptu_prediction.tsv"), emit: ptu
    tuple val(meta), val(plasmid_name), path("$plasmid_name/copla.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${plasmid_name}

    copla \\
         "$file" \\
         "${copladb}/Copla_RS84/RS84f_sHSBM.pickle" \\
         "${copladb}/Copla_RS84/CoplaDB.fofn" \\
         ${plasmid_name} > ${plasmid_name}/copla.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$( copla --version | cut -f2 -d' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${plasmid_name}
    touch ${plasmid_name}/copla.txt
    touch ${plasmid_name}/stub.ptu_prediction.tsv
    touch ${plasmid_name}/stub.qry_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copla: \$( copla --version | cut -f2 -d' ')
    END_VERSIONS
    """
}
