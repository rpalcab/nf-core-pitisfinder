/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'

// include { SAMPLESUMMARY          } from '../modules/local/samplesummary/main'
include { MERGE_ANNOTATIONS      } from '../modules/local/mergeannotations/main'
include { ISESCAN                } from '../modules/local/isescan/main'

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

    RVD_ANNOTATION.out.amr_report
    RVD_ANNOTATION.out.vf_report

    // MERGE GBK AND AMR ANNOTATIONS
    ch_samplesheet.map {  meta, fasta, gbk ->
        return [ meta, gbk ]
        }
        .join ( RVD_ANNOTATION.out.amr_report )
        .join ( RVD_ANNOTATION.out.vf_report )
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
                            summary ? [ meta, summary_list + [ summary ], gbk_list ] : [ meta, file_list, summary_list, gbk_list ]
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
              .set { ch_samplesummary }

    SAMPLE_SUMMARY(ch_samplesummary)

    // SAMPLESUMMARY(ch_samplesummary)


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
