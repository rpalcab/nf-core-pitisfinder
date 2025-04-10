process PHIGARO_SETUP {
    tag "phigaro_single"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phigaro:2.4.0--pyhdfd78af_0':
        'biocontainers/phigaro:2.4.0--pyhdfd78af_0' }"

    input:
    path(pvog_input)

    output:
    path("config.yaml") , emit: config
    path("pvog")        , emit: pvog
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // WARNING: Ensure to update database version if it is updated.
    def pvog_version = "18-Aug-2017"
    """
    set -e  # Exit immediately if a command exits with a non-zero status

    if [ -z "${pvog_input}" ]; then
        mkdir pvog
        # Function to download file and check if it was successful
        download_file() {
            local file=\$1
            local max_attempts=3
            local attempt=1

            while [ \$attempt -le \$max_attempts ]; do
                if wget -q http://download.ripcm.com/phigaro/\$file -O pvog/\$file; then
                    echo "Successfully downloaded \$file"
                    return 0
                else
                    echo "Attempt \$attempt failed to download \$file"
                    attempt=\$((attempt + 1))
                    sleep 5  # Wait for 5 seconds before retrying
                fi
            done

            echo "Failed to download \$file after \$max_attempts attempts"
            return 1
        }

        # Download files with error handling
        for file in allpvoghmms allpvoghmms.h3f allpvoghmms.h3i allpvoghmms.h3m allpvoghmms.h3p; do
            if ! download_file \$file; then
                echo "Error: Failed to download \$file. Exiting."
                exit 1
            fi
        done
    else
        # Check if the input directory is already named 'pvog'
        if [ "\$(basename ${pvog_input})" != "pvog" ]; then
            mv ${pvog_input} pvog
        fi
    fi

    cat <<EOF > config.yaml
    hmmer:
     bin: hmmsearch
     e_value_threshold: 0.00445
     pvog_path: ./pvog/allpvoghmms
    phigaro:
     mean_gc: 0.46354823199323625
     penalty_black: 2.2
     penalty_white: 0.7
     threshold_max_abs: 52.96
     threshold_max_basic: 46.0
     threshold_max_without_gc: 11.42
     threshold_min_abs: 50.32
     threshold_min_basic: 45.39
     threshold_min_without_gc: 11.28
     window_len: 32
    prodigal:
     bin: prodigal
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phigaro_db: \$(echo "${pvog_version}")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    // WARNING: Ensure to update database version if it is updated.
    def pvog_version = "18-Aug-2017"
    """
    mkdir pvog
    touch config.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phigaro_db: \$(echo "${pvog_version}")
    END_VERSIONS
    """
}