/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'
include { SAMPLESUMMARY          } from '../modules/local/samplesummary/main'
include { MERGE_ANNOTATIONS      } from '../modules/local/mergeannotations/main'
include { ISESCAN                } from '../modules/local/isescan/main'

include { RVD_ANNOTATION         } from '../subworkflows/local/rvd_annotation'
include { PLASMID_ANALYSIS       } from '../subworkflows/local/plasmid_analysis'
include { INTEGRON_ANALYSIS      } from '../subworkflows/local/integron_analysis'
include { PROPHAGE_ANALYSIS      } from '../subworkflows/local/prophage_analysis'
include { ICE_ANALYSIS           } from '../subworkflows/local/ice_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    //
    // INITIALIZE CHANNELS
    //
    ch_versions = Channel.empty()

    // INITIALIZE SUMMARY CHANNEL
    ch_samplesheet.map { meta, fasta, faa, gbk ->
        return [ meta, [] ]
        }.set { ch_summary }

    // PREPARE ONLY FASTA CHANNEL
    ch_samplesheet.map { meta, fasta, faa, gbk ->
        return [ meta, fasta ]
        }.set { ch_fasta }

    // PREPARE FAA CHANNEL
    ch_samplesheet
        .map { meta, fasta, faa, gbk ->
            return [ meta, faa ]
        }.set { ch_faa }

    // RES, VIR, DEF ANNOTATION
    RVD_ANNOTATION (
            ch_fasta,
            params.df_db ? params.df_db : null
        )

    RVD_ANNOTATION.out.amr_report
    RVD_ANNOTATION.out.vf_report

    // MERGE GBK AND AMR ANNOTATIONS
    ch_samplesheet.map {  meta, fasta, faa, gbk ->
        return [ meta, gbk ]
        }
        .join ( RVD_ANNOTATION.out.amr_report )
        .join ( RVD_ANNOTATION.out.vf_report )
        .set { ch_mergeann }
    MERGE_ANNOTATIONS ( ch_mergeann )

    // PREPARE FULL CHANNEL WITH MERGED ANNOTATIONS
    ch_samplesheet
        .join(MERGE_ANNOTATIONS.out.gbk)
        .map {  meta, fasta, faa, gbk, merged_gbk ->
            return [ meta, fasta, faa, merged_gbk ]
        }.set { ch_full }

    // PREPARE ONLY GBK WITH MERGED ANNOTATIONS CHANNEL
     ch_full
        .map { meta, fasta, faa, gbk ->
            return [ meta, gbk ]
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
                        .map { meta, file_list, summary ->
                            summary ? [ meta, file_list + [ summary ] ] : [ meta, file_list ]
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
                        .map { meta, file_list, summary ->
                            summary ? [ meta, file_list + [ summary ] ] : [ meta, file_list ]
                    }
    }

    //
    // IS
    //
    if ( !params.skip_is ) {
        ISESCAN (ch_fasta)
        ch_versions = ch_versions.mix( ISESCAN.out.versions )

        ch_summary = ch_summary
                        .join(ISESCAN.out.summary, remainder: true)
                        .map { meta, file_list, summary ->
                            summary ? [ meta, file_list + [summary] ] : [ meta, file_list ]
                    }
    }

    //
    // PROPHAGES
    //
    if ( !params.skip_prophages ) {
        PROPHAGE_ANALYSIS (
                    ch_gbk,
                    params.phispy_db ? params.phispy_db : null
                )
        ch_versions = ch_versions.mix(PROPHAGE_ANALYSIS.out.versions)

        ch_summary = ch_summary
                        .join(PROPHAGE_ANALYSIS.out.summary, remainder: true)
                        .map { meta, file_list, summary ->
                            summary ? [ meta, file_list + [ summary ] ] : [ meta, file_list ]
                    }
    }

    //
    // ICEs
    //
    if ( !params.skip_ices ) {
        ICE_ANALYSIS(ch_faa, ch_gbk, PLASMID_ANALYSIS.out.contig_report)
        ch_versions = ch_versions.mix(ICE_ANALYSIS.out.versions)
    }

    //
    // SAMPLE SUMMARY
    //
    ch_summary.join( MERGE_ANNOTATIONS.out.gbk )
              .set { ch_samplesummary }

    SAMPLESUMMARY(ch_samplesummary)


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
