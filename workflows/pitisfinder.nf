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
include { IS_BLAST                    } from '../modules/local/isblast/main'
include { IS_PARSER                   } from '../modules/local/isparser/main'
include { ISESCAN                     } from '../modules/local/isescan/main'
include { PHISPY                      } from '../modules/nf-core/phispy/main'
include { MACSYFINDER                 } from '../modules/local/macsyfinder/macsyfinder/main'
include { VERIFYMODEL                 } from '../modules/local/macsyfinder/verifymodel/main'
include { ICEFINDER2                  } from '../modules/local/icefinder2/main'
include { SAMPLESUMMARY               } from '../modules/local/samplesummary/main'

include { PLASMID_ANALYSIS            } from '../subworkflows/local/plasmid_analysis'
include { INTEGRON_ANALYSIS           } from '../subworkflows/local/integron_analysis'

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

    if ( !params.skip_is ) {
        ////
        //// INSERTION SEQUENCES (ESTO IRÁ A SUBWORKFLOW)
        ////
        // ch_is_input = Channel.empty()
        // if ( !params.skip_plasmids ) {
        //     ch_is_input = ch_mobsuite_chr
        // } else {
        //     ch_is_input = ch_fasta
        // }
        // ch_isdb = Channel.empty()
        // if (params.is_db){
        //     ch_isdb = Channel.value(params.is_db)
        //     // BLASTn search
        //     IS_BLAST (ch_is_input, ch_isdb)
        //     ch_versions = ch_versions.mix( IS_BLAST.out.versions )
        //     // Results filtering
        //     ch_raw_is = IS_BLAST.out.report
        //     IS_PARSER (ch_raw_is)
        // }
        ISESCAN (ch_fasta)
        ch_versions = ch_versions.mix( ISESCAN.out.versions )
        ch_summary = ch_summary
                        .join(ISESCAN.out.summary, remainder: true)
                        .map { meta, file_list, summary ->
                            summary ? [meta, file_list + [summary]] : [meta, file_list]
                    }
    }

    if ( !params.skip_prophages ) {
        //
        // PHIPSY
        //
        ch_phispydb = Channel.value([])
        if (params.phispy_db){
            ch_phispydb = Channel.value(params.phispy_db)
        }
        ch_full.map { meta, fasta, faa, gbk ->
            return [ meta, gbk ]
            }.set { ch_phispy }
        PHISPY (ch_phispy, ch_phispydb)
        ch_versions = ch_versions.mix( PHISPY.out.versions )

        PROCESS_PHISPY (PHISPY.out.phage_gbk, PHISPY.out.phage_fasta, PHISPY.out.prophage_tsv)
        PROCESS_PHISPY.out.gbk
            .flatMap { meta, gbk_files ->
            if (gbk_files instanceof List) {
                // If gbk_files is a list, process each file
                return gbk_files.collect { gbk_file ->
                    def pp_meta = [id: gbk_file.name.tokenize('.')[0]]
                    [meta, pp_meta, gbk_file]
                }
            } else {
                // If gbk_files is a single file, process it directly
                def pp_meta = [id: gbk_files.name.tokenize('.')[0]]
                return [[meta, pp_meta, gbk_files]]
            }
        }.set { ch_pplin }
        ch_summary = ch_summary
                        .join(PROCESS_PHISPY.out.summary, remainder: true)
                        .map { meta, file_list, summary ->
                            summary ? [meta, file_list + [summary]] : [meta, file_list]
                    }
    }

    if ( !params.skip_ices ) {
        //
        // ICEs (ESTO IRÁ A SUBWORKFLOW)
        //
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
