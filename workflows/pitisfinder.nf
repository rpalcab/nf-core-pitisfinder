/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'
include { VISUALIZE_CIRCULAR          } from '../modules/local/visualize/circular/main'
include { VISUALIZE_LINEAR            } from '../modules/local/visualize/linear/main'
include { ISESCAN                     } from '../modules/local/isescan/main'
include { PHISPY                      } from '../modules/nf-core/phispy/main'
include { MACSYFINDER                 } from '../modules/local/macsyfinder/macsyfinder/main'
include { VERIFYMODEL                 } from '../modules/local/macsyfinder/verifymodel/main'
include { ICEFINDER2                  } from '../modules/local/icefinder2/main'
include { SAMPLESUMMARY               } from '../modules/local/samplesummary/main'

include { PLASMID_ANALYSIS            } from '../subworkflows/local/plasmid_analysis'
include { INTEGRON_ANALYSIS           } from '../subworkflows/local/integron_analysis'
include { PROPHAGE_ANALYSIS           } from '../subworkflows/local/prophage_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_FAA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(faa), path(gbk), path(chr)

    output:
    tuple val(meta), path("${meta.id}/*_filtered.faa"), emit: chr

    script:
    def prefix = "${meta.id}"
    """
    cut -f2,5 $chr | grep "chromosome" | cut -f2 > chr.txt
    faa_contig_filter.py -f $faa -g $gbk -c chr.txt -o $prefix
    """
}

process MERGE_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(amr), path(gbk)

    output:
    tuple val(meta), path("${meta.id}_merged.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    """
    merge_amr.py -t $amr -g $gbk -o "$prefix"_merged.gbk
    """
}

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    // INITIALIZE SUMMARY CHANNEL
    ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        return [ meta, [] ]
        }.set { ch_summary }

    // PREPARE ONLY FASTA CHANNEL
    ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        return [ meta, fasta ]
        }.set { ch_fasta }

    // MERGE GBK AND AMR ANNOTATIONS
    ch_samplesheet.map {  meta, fasta, faa, gbk, amr ->
        return [ meta, amr, gbk ]
        }.set { ch_mergeann }
    MERGE_ANNOTATIONS ( ch_mergeann )

    // PREPARE FULL CHANNEL WITH MERGED ANNOTATIONS
    ch_samplesheet
        .join(MERGE_ANNOTATIONS.out.gbk)
        .map {  meta, fasta, faa, gbk, amr, merged_gbk ->
            return [ meta, fasta, faa, merged_gbk ]
        }.set { ch_full }

    // PREPARE ONLY GBK WITH MERGED ANNOTATIONS CHANNEL
     ch_full
        .map { meta, fasta, faa, gbk ->
            return [ meta, gbk ]
        }.set { ch_gbk }

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

        // MACSYFINDER
        ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
            return [ meta, faa, gbk ]
        }.join(PLASMID_ANALYSIS.out.contig_report)
        .set { ch_msy_preproc }
        FILTER_FAA (ch_msy_preproc)
        ch_msymodel = Channel.value('CONJScan/Plasmids')
        VERIFYMODEL (ch_msymodel)
        ch_model = VERIFYMODEL.out.model
        MACSYFINDER (FILTER_FAA.out.chr, ch_model)
        ch_versions = ch_versions.mix( MACSYFINDER.out.versions )

        ////
        //// ICEFINDER2
        ////
        // ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        //     return [ meta, gbk ]
        // }.set { ch_icefinder }
        // ICEFINDER2 ( ch_icefinder )
        // ch_versions = ch_versions.mix( ICEFINDER2.out.versions )
    }

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
