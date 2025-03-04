process COPLA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/copla:latest' :
        'docker.io/rpalcab/copla:latest' }"

    input:
    tuple val(meta), path(plasmid_files)

    output:
    path "copla_results/${meta.id}/", emit: copla_results
    path "copla.txt", emit: copla_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    # Descargar bases de datos si no existen
    if [ ! -d "databases" ]; then
        bin/download_Copla_databases.sh
    fi

    mkdir -p copla_results/$prefix

    # Filtrar plásmidos <600kb y >1kb
    for i in \$(find ${plasmid_files} -size -600k -size +1k); do
        new_name=\$(basename "\$i" | rev | cut -f1 -d'/' | cut -f2- -d'.' | cut -f1 -d'_' | rev)
        cp "\$i" 05_plasmids/\${new_name}_$prefix.fasta

        echo "Sample: $prefix" >> copla.txt
        echo "Contig: \${i:(-11):5}" >> copla.txt

        python3 bin/copla.py \\
            "\$i" \\
            databases/Copla_RS84/RS84f_sHSBM.pickle \\
            databases/Copla_RS84/CoplaDB.fofn \\
            copla_results/$prefix
    done

    python3 bin/copla.py ${i} databases/Copla_RS84/RS84f_sHSBM.pickle databases/Copla_RS84/CoplaDB.fofn 08_Anotacion/${j}/copla

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
