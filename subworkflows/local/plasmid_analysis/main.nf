include { MOBSUITE_RECON      } from '../../../modules/nf-core/mobsuite/recon/main'
include { COPLA_COPLADBDOWNLOAD     } from '../../../modules/local/copla/copladbdownload/main'
include { COPLA_COPLA     } from '../../../modules/local/copla/copla/main'
include { PLASMID_PARSER     } from '../../../modules/local/plasmidparser/main'
include { VISUALIZE_CIRCULAR     } from '../../../modules/local/visualize/circular/main'

// TODO: Move these processes to plasmid_utils or slt
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

process PLASMID_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("plasmid_summary.tsv"), emit: summary

    script:
    def prefix = "${meta.id}"
    """
    echo -e "Contig\tName\tStart\tEnd\tAMR" > plasmid_summary.tsv
    tail -q -n+2 $report | cut -f2,6,3,4,19 | awk -F'\t' 'BEGIN{OFS="\t"} {print \$1, \$2":"\$3, 1, \$4, \$5}' >> plasmid_summary.tsv
    """
}

workflow PLASMID_ANALYSIS {

    take:
    ch_fasta        // channel (mandatory): [ val(meta), fasta ]
    ch_gbk          // channel(mandatory): [ val(meta), gbk]
    copla_db        // path (optional): copla_db

    main:

    ch_versions = Channel.empty()

    // MOBSUITE
    MOBSUITE_RECON ( ch_fasta )
    ch_versions = ch_versions.mix( MOBSUITE_RECON.out.versions )

    // RENAME PLASMIDS
    ch_rename = RENAME_PLASMIDS ( MOBSUITE_RECON.out.plasmids )
    ch_rename
            .flatMap { meta, plasmid_files ->
                def file_list = plasmid_files instanceof List ? plasmid_files : [ plasmid_files ]
                file_list.collect { file ->
                    def plasmid_name = file.baseName
                    return tuple( meta, plasmid_name, file )
                }
            }
            .set{ ch_plasmids }

    // COPLA
    ch_copladb = Channel.empty()
    if (!copla_db){
        COPLA_COPLADBDOWNLOAD ()
        ch_copladb = COPLA_COPLADBDOWNLOAD.out.db
    } else {
        ch_copladb = Channel.value(file(copla_db))
    }
    COPLA_COPLA ( ch_plasmids, ch_copladb )
    ch_versions = ch_versions.mix( COPLA_COPLA.out.versions )

    // PLASMID_PARSER
    COPLA_COPLA.out.query
        .join(COPLA_COPLA.out.ptu, by: [ 0, 1 ])
        .set { ch_coplajoint }

    ch_gbk
        .join(MOBSUITE_RECON.out.mobtyper_results)
        .join(MOBSUITE_RECON.out.contig_report)
        .set { ch_mobsample }

    ch_mobsample
        .cross( ch_coplajoint )
        .map { mob, copla ->
            def meta = mob[0]
            def gbk = mob[1]
            def mob_typer = mob[2]
            def contig_report = mob[3]
            def plasmid_name = copla[1]
            def qry = copla[2]
            def ptu = copla[3]
            return [ meta, plasmid_name, qry, ptu, gbk, mob_typer, contig_report ]
        }
        .set { ch_plasmidparser }

    PLASMID_PARSER( ch_plasmidparser )

    // CREATE SUMMARY
    PLASMID_PARSER.out.report
        .map { meta, plasmid, report -> [ meta, report ] }
        .groupTuple()
        .map { meta, report -> [ meta, report.flatten() ] }
        .set { ch_plasmidsummary }

    PLASMID_SUMMARY( ch_plasmidsummary )
    VISUALIZE_CIRCULAR( PLASMID_PARSER.out.gbk )

    emit:
    plasmids        = MOBSUITE_RECON.out.plasmids       // channel: [ val(meta), [ plasmid_fastas ] ]
    chromosome      = MOBSUITE_RECON.out.chromosome     // channel: [ val(meta), [ chromosome_fastas ] ]
    contig_report   = MOBSUITE_RECON.out.contig_report  // channel: [ val(meta), contig_report ]
    summary         = PLASMID_SUMMARY.out.summary       // channel: [ val(meta), summary ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}

