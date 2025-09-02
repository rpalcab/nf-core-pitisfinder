/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'

include { MERGE_ANNOTATIONS      } from '../modules/local/mergeannotations/main'
include { ISESCAN                } from '../modules/local/isescan/main'
include { VISUALIZE_PLASMID      } from '../modules/local/visualize/plasmid/main'
include { MGESUMMARY            } from '../modules/local/mgesummary/main'

include { RVD_ANNOTATION         } from '../subworkflows/local/rvd_annotation'
include { PLASMID_ANALYSIS       } from '../subworkflows/local/plasmid_analysis'
include { INTEGRON_ANALYSIS      } from '../subworkflows/local/integron_analysis'
include { PROPHAGE_ANALYSIS      } from '../subworkflows/local/prophage_analysis'
include { ICE_ANALYSIS           } from '../subworkflows/local/ice_analysis'
include { SAMPLE_SUMMARY         } from '../subworkflows/local/summary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INITIALIZE CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = Channel.empty()

    // INITIALIZE SUMMARY CHANNEL
    ch_samplesheet.map { meta, fasta, gbk ->
        return [ meta, [], [] ]
        }.set { ch_summary }

    // PREPARE ONLY FASTA CHANNEL
    ch_samplesheet.map { meta, fasta, gbk ->
        return [ meta, fasta ]
        }.set { ch_fasta }


    // RES, VIR, DEF ANNOTATION
    RVD_ANNOTATION (
            ch_fasta,
            params.df_db ? params.df_db : null
        )

    // MERGE GBK AND AMR ANNOTATIONS
    ch_samplesheet.map {  meta, fasta, gbk ->
        return [ meta, gbk ]
        }
        .join ( RVD_ANNOTATION.out.amr_report )
        .join ( RVD_ANNOTATION.out.vf_report )
        .join ( RVD_ANNOTATION.out.df_report )
        .set { ch_mergeann }
    MERGE_ANNOTATIONS ( ch_mergeann )

    // PREPARE CHANNEL WITH MERGED ANNOTATIONS
    ch_samplesheet
        .join(MERGE_ANNOTATIONS.out.gbk)
        .map {  meta, fasta, gbk, merged_gbk ->
            return [ meta, merged_gbk ]
        }.set { ch_gbk }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START EGM SECTION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    //
    // PLASMIDS
    //
    ch_plasmid_contig_report = Channel.empty()

    if ( !params.skip_plasmids ) {
        PLASMID_ANALYSIS(
                    ch_fasta,
                    ch_gbk,
                    params.copla_db ? params.copla_db : null
                )
        ch_versions = ch_versions.mix( PLASMID_ANALYSIS.out.versions )

        ch_summary = ch_summary
                        .join( PLASMID_ANALYSIS.out.summary, remainder: true )
                        .join( PLASMID_ANALYSIS.out.genomic_gbk, remainder: true )
                        .map { meta, summary_list, gbk_list, summary, genomic_gbk ->
                            summary ? [ meta, summary_list + [ summary ], gbk_list + [ genomic_gbk ] ] : [ meta, summary_list, gbk_list ]
                    }

        ch_plasmid_contig_report = PLASMID_ANALYSIS.out.contig_report
    } else {
        // Create empty channel with same structure when plasmids are skipped
        ch_plasmid_contig_report = ch_samplesheet.map { meta, fasta, gbk ->
            return [ meta, "" ]
        }
    }

    //
    // INTEGRONS
    //
    if ( !params.skip_integrons ) {
        INTEGRON_ANALYSIS(
                    ch_fasta,
                    ch_gbk
                )
        ch_versions = ch_versions.mix(INTEGRON_ANALYSIS.out.versions)

        ch_summary = ch_summary
                        .join( INTEGRON_ANALYSIS.out.summary, remainder: true )
                        .join( INTEGRON_ANALYSIS.out.genomic_gbk, remainder: true )
                        .map { meta, summary_list, gbk_list, summary, genomic_gbk ->
                            summary ? [ meta, summary_list + [ summary ], gbk_list + [ genomic_gbk ] ] : [ meta, summary_list, gbk_list ]
                    }
    }

    //
    // IS
    //
    if ( !params.skip_is ) {
        ISESCAN (ch_fasta)
        ch_versions = ch_versions.mix( ISESCAN.out.versions )

        ch_summary = ch_summary
                        .join( ISESCAN.out.summary, remainder: true )
                        .map { meta, summary_list, gbk_list, summary ->
                            summary ? [ meta, summary_list + [ summary ], gbk_list ] : [ meta, summary_list, gbk_list ]
                    }
    }

    //
    // PROPHAGES
    //
    if ( !params.skip_prophages ) {
        PROPHAGE_ANALYSIS (
                    ch_fasta,
                    ch_gbk,
                    params.genomad_db ? params.genomad_db : null
                )
        ch_versions = ch_versions.mix(PROPHAGE_ANALYSIS.out.versions)

        ch_summary = ch_summary
                        .join( PROPHAGE_ANALYSIS.out.summary, remainder: true )
                        .join( PROPHAGE_ANALYSIS.out.genomic_gbk, remainder: true )
                        .map { meta, summary_list, gbk_list, summary, genomic_gbk ->
                            summary ? [ meta, summary_list + [ summary ], gbk_list + [ genomic_gbk ] ] : [ meta, summary_list, gbk_list ]
                    }
    }

    //
    // ICEs
    //
    // if ( !params.skip_ices ) {
    //     ICE_ANALYSIS(ch_faa, ch_gbk, PLASMID_ANALYSIS.out.contig_report)
    //     ch_versions = ch_versions.mix(ICE_ANALYSIS.out.versions)
    // }

    //
    // SAMPLE SUMMARY
    //
    ch_summary.join( ch_gbk )
              .join( ch_plasmid_contig_report )
              .set { ch_samplesummary }

    SAMPLE_SUMMARY(ch_samplesummary)

    //
    // PLASMID VISUALIZATION
    //
    if ( !params.skip_plasmids ) {
        PLASMID_ANALYSIS.out.plasmid_gbk
            .groupTuple(by: 0)
            .join( SAMPLE_SUMMARY.out.tsv )
            .transpose()
            .map { meta, plasmid_name, gbk_file, summary_tsv ->
                return [meta, plasmid_name, gbk_file, summary_tsv]
            }
            .set { ch_visual }

        VISUALIZE_PLASMID( ch_visual )
    }

    //
    // FINAL MGE SUMMARY
    //
    // Gral
    ch_mgesum_gral_tsv = SAMPLE_SUMMARY.out.tsv
                            .map { meta, tsv ->
                                return [ tsv ]
                            }
                            .collect()
    ch_mgesum_gral_png = SAMPLE_SUMMARY.out.png
                            .map { meta, name, png ->
                                return [ png ]
                            }
                            .collect()

    // Plasmids
    ch_mgesum_plasmids_png = params.skip_plasmids ? Channel.empty() : VISUALIZE_PLASMID.out.png
                                                .map { meta, pl_id, png ->
                                                    return [ png ]
                                                }
                                                .collect()

    ch_mgesum_plasmids_reports = params.skip_plasmids ? Channel.empty() : PLASMID_ANALYSIS.out.plasmid_report
                                                .map { meta, name, tsv ->
                                                    return [ tsv ]
                                                }
                                                .collect()

    // Integrons
    ch_mgesum_integrons_tsv = params.skip_integrons ? Channel.empty() : INTEGRON_ANALYSIS.out.summary
                                                .map { meta, tsv ->
                                                    return [ tsv ]
                                                }
                                                .collect()

    ch_mgesum_integrons_png = params.skip_integrons ? Channel.empty() : INTEGRON_ANALYSIS.out.png
                                                .map { meta, int_id, png ->
                                                    return [ png ]
                                                }
                                                .collect()


    // Prophages
    ch_mgesum_prophages_tsv = params.skip_prophages ? Channel.empty() : PROPHAGE_ANALYSIS.out.summary
                                                .map { meta, tsv ->
                                                    return [ tsv ]
                                                }
                                                .collect()

    ch_mgesum_prophages_png = params.skip_prophages ? Channel.empty() : PROPHAGE_ANALYSIS.out.png
                                                .map { meta, prophage_id, png ->
                                                    return [ png ]
                                                }
                                                .collect()

    // ch_mgesum_gral_tsv.view()
    // ch_mgesum_gral_png.view()
    // ch_mgesum_plasmids_tsv.view()
    // ch_mgesum_plasmids_png.view()
    // ch_mgesum_integrons_tsv.view()
    // ch_mgesum_integrons_png.view()
    // ch_mgesum_prophages_tsv.view()
    // ch_mgesum_prophages_png.view()

    MGESUMMARY(
            ch_mgesum_gral_tsv.ifEmpty([]),
            ch_mgesum_gral_png.ifEmpty([]),
            ch_mgesum_plasmids_png.ifEmpty([]),
            ch_mgesum_plasmids_reports.ifEmpty([]),
            ch_mgesum_integrons_tsv.ifEmpty([]),
            ch_mgesum_integrons_png.ifEmpty([]),
            ch_mgesum_prophages_tsv.ifEmpty([]),
            ch_mgesum_prophages_png.ifEmpty([])
            )
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'pitisfinder_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
