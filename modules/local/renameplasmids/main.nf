process RENAME_PLASMIDS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(plasmid_files)

    output:
    tuple val(meta), path("*.fasta"), emit: filtered_plasmids, optional: true

    script:
    """
    for i in ${plasmid_files}; do
        new_name=\$(basename "\$i" | sed 's/plasmid_//' | sed 's/.fasta//')
        mv "\$i" \${new_name}_${meta.id}.fasta
    done
    """
}
