process ICEBERG_DB_DOWNLOAD {
    tag "iceberg_dwnld"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--h66d330f_4':
        'biocontainers/blast:2.16.0--h66d330f_4' }"

    output:
    path("ICE_db.fasta")                    , emit: db
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    wget https://tool2-mml.sjtu.edu.cn/ICEberg3/data/download/ICE_seq_all.fas
    # Check all characters in headers are valid for makeblastdb
    sed -i 's/ .*//g' ICE_seq_all.fas
    awk '{
    if (\$0 ~ /^>/) {
        header = substr(\$0, 2)
        gsub(/[^A-Za-z0-9_.-]/, "_", header)
        print ">" header
    } else {
        print
        }
    }' ICE_seq_all.fas > ICE_db.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | grep "GNU Wget" | sed 's/GNU Wget //; s/built.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ICE_db.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | grep "GNU Wget" | sed 's/GNU Wget //; s/built.*\$//')
    END_VERSIONS
    """
}
