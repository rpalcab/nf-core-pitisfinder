include { MOBSUITE_RECON        } from '../../../modules/nf-core/mobsuite/recon/main'
include { PLASMIDMARKERS        } from '../../../modules/local/plasmidmarkers/main'
include { COPLA_COPLADBDOWNLOAD } from '../../../modules/local/copla/copladbdownload/main'
include { COPLA_COPLA           } from '../../../modules/local/copla/copla/main'
include { PLASMID_PARSER        } from '../../../modules/local/plasmidparser/main'
include { RENAME_PLASMIDS       } from '../../../modules/local/renameplasmids/main'
include { PLASMID_SUMMARY       } from '../../../modules/local/plasmidsummary/main'
include { VISUALIZE_PLASMID     } from '../../../modules/local/visualize/plasmid/main'

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

    // UPDATE GENOMIC GBK
    MOBSUITE_RECON.out.biomarkers
                .join(ch_gbk)
                .set { ch_updategbk }

    // ANNOTATE PLASMID MARKERS
    PLASMIDMARKERS(ch_updategbk)

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

    PLASMIDMARKERS.out.gbk
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

    // // VISUALIZATION
    // PLASMID_PARSER.out.gbk.map
    //             { meta, plasmid_name, gbk ->
    //             return [ meta, plasmid_name, gbk ]
    //             }.set { ch_visual }
    // VISUALIZE_PLASMID( ch_visual )

    // CREATE SUMMARY
    PLASMID_PARSER.out.report
            .groupTuple()
            .map { meta, plasmid, report ->
                return [ meta, report.flatten() ]
            }
            .set { ch_reports }

    PLASMID_PARSER.out.gbk
            .groupTuple()
            .map { meta, plasmid, gbk ->
                return [ meta, gbk.flatten() ]
            }
            .set { ch_plasmidgbks }

    ch_reports
            .join(ch_plasmidgbks)
            .join(PLASMIDMARKERS.out.gbk)
            .set { ch_plasmidsummary }

    PLASMID_SUMMARY( ch_plasmidsummary )


    emit:
    plasmids        = MOBSUITE_RECON.out.plasmids       // channel: [ val(meta), [ plasmid_fastas ] ]
    chromosome      = MOBSUITE_RECON.out.chromosome     // channel: [ val(meta), [ chromosome_fastas ] ]
    contig_report   = MOBSUITE_RECON.out.contig_report  // channel: [ val(meta), contig_report ]
    summary         = PLASMID_SUMMARY.out.summary       // channel: [ val(meta), summary ]
    genomic_gbk     = PLASMIDMARKERS.out.gbk            // channel: [ val(meta), gbk ]
    plasmid_gbk     = PLASMID_PARSER.out.gbk            // channel: [ val(meta), plasmid_name, gbk ]
    versions        = ch_versions                       // channel: [ versions.yml ]
}

