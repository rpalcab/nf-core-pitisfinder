process PVOGDOWNLOAD {
    tag "pvog_download"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4':
        'biocontainers/blast:2.16.0--h66d330f_4' }"

    output:
    path "db/pVOGs.hmm"           , emit: pvogs_db
    path "pvogs_annotations.csv"  , emit: pvogs_ann

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def hmm_profile = "https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"
    def pvog_annotations = "https://raw.githubusercontent.com/bobeobibo/phigaro/master/phigaro/to_html/pvogs_annotations.csv"
    """
    download_file() {
        local file=\$1
        local max_attempts=3
        local attempt=1

        while [ \$attempt -le \$max_attempts ]; do
            if wget -q \$file; then
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

    for file in $hmm_profile $pvog_annotations ; do
        if ! download_file \$file; then
            echo "Error: Failed to download \$file. Exiting."
            exit 1
        fi
    done

    tar -xzf AllvogHMMprofiles.tar.gz && rm AllvogHMMprofiles.tar.gz
    cat AllvogHMMprofiles/* > pVOGs.hmm
    rm -r AllvogHMMprofiles/
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir pvog
    touch pvog/pVOGs.hmm
    touch pvogs_annotations.csv
    """
}
