process MGESUMMARY {
    tag 'final_summary'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/jinja2_numpy_pandas_plotly_python:ae2b2f0af78b69ae':
        'community.wave.seqera.io/library/jinja2_numpy_pandas_plotly_python:ae2b2f0af78b69ae' }"

    input:
    path(gral_tsv)
    path(gral_png)
    path(plasmids_png)
    path(plasmids_reports)
    path(integrons_tsv)
    path(integrons_png)
    path(prophages_tsv)
    path(prophages_png)

    output:
    path("mge_summary.html"), emit: html
    path("images")          , emit: images

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gral_tsv_str = gral_tsv.join(',')
    def gral_tsv_opt = gral_tsv_str ? "-g ${gral_tsv_str}" : ""
    def gral_png_str = gral_png.join(',')
    def gral_png_opt = gral_png_str ? "-G ${gral_png_str}" : ""
    def plasmids_png_str = plasmids_png.join(',')
    def plasmids_png_opt = plasmids_png_str ? "-P ${plasmids_png_str}" : ""
    def plasmids_reports_str = plasmids_reports.join(',')
    def plasmids_reports_opt = plasmids_reports_str ? "--plasmid_reports ${plasmids_reports_str}" : ""
    def integrons_tsv_str = integrons_tsv.join(',')
    def integrons_tsv_opt = integrons_tsv_str ? "-i ${integrons_tsv_str}" : ""
    def integrons_png_str = integrons_png.join(',')
    def integrons_png_opt = integrons_png_str ? "-I ${integrons_png_str}" : ""
    def prophages_tsv_str = prophages_tsv.join(',')
    def prophages_tsv_opt = prophages_tsv_str ? "-f ${prophages_tsv_str}" : ""
    def prophages_png_str = prophages_png.join(',')
    def prophages_png_opt = prophages_png_str ? "-F ${prophages_png_str}" : ""
    """
    final_report.py \\
        ${gral_tsv_opt} \\
        ${gral_png_opt} \\
        ${plasmids_reports_opt} \\
        ${plasmids_png_opt} \\
        ${integrons_tsv_opt} \\
        ${integrons_png_opt} \\
        ${prophages_tsv_opt} \\
        ${prophages_png_opt} \\
        -o mge_summary.html
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch mge_summary.html

    final_report.py \\
        -i ${results_dir} \\
        -o mge_summary.html

    """
}
