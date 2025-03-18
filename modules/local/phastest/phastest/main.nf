process PHASTEST_PHASTEST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/phastest_single:1.1' :
        'docker.io/rpalcab/phastest_single:1.1' }"

    input:
    tuple val(meta), path(fasta)
    path (phastestdb)

    output:
    tuple val(meta), path("${meta.id}"), emit: outdir
    // No output version available

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    """
    phastest -i fasta -s "$fasta" -d "$phastestdb" --phage-only --yes 
    """

    stub:
    """
    mkdir -p ${prefix}
    """
}
