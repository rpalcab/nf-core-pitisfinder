process PROCESS_PHISPY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://rpalcab/visualizer:1.0':
        'docker.io/rpalcab/visualizer:1.0' }"

    input:
    tuple val(meta), path(gbk)
    tuple val(meta), path(fasta)
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("pp_*_${meta.id}.gbk"), emit: gbk
    tuple val(meta), path("pp_*_${meta.id}.fasta"), emit: fasta
    tuple val(meta), path("pp_*_${meta.id}.png"), emit: png
    tuple val(meta), path("prophage_summary.tsv"), emit: summary

    script:
    def prefix = "${meta.id}"
    """
    # Split gbk into individual prophages
    csplit -z $gbk /^LOCUS/ '{*}' -f tmp_pp_ -b %d_${prefix}.gbk
    # Rename gbk files starting from 1
    n=1
    for file in \$(ls tmp_pp_*_${prefix}.gbk | sort -V); do
        newname=\$(printf "pp_%d_${prefix}.gbk" "\$n")
        mv "\$file" "\$newname"
        linear_plot.py -i \$newname -o \${newname%.*}.png
        n=\$((n + 1))
    done

    # Split fasta into individual prophages
    csplit -z $fasta /^\\>/ '{*}' -f tmp_pp_ -b %d_${prefix}.fasta
    # Rename fasta files starting from 1
    n=1
    for file in \$(ls tmp_pp_*_${prefix}.fasta | sort -V); do
        newname=\$(printf "pp_%d_${prefix}.fasta" "\$n")
        mv "\$file" "\$newname"
        n=\$((n + 1))
    done

    # Format summary table
    awk -v s="$prefix" 'NR==1 {print \$0; next} { \$1 = \$1 "_" s; print }' OFS='\\t' $tsv > prophage_summary.tsv
    sed -i -e 's/Prophage number/Name/' -e 's/Stop/End/' prophage_summary.tsv
    """
}