include { MERGEEGM    } from '../../../modules/local/mergeegm/main'
include { MERGEGBKS   } from '../../../modules/local/mergegbks/main'
include { VISUALIZE_CIRCULAR    } from '../../../modules/local/visualize/circular/main'

workflow SAMPLE_SUMMARY {
    take:
    ch_input    // channel: [meta, tsv_list, gbk_list, gbk]

    main:
    // ch_input.view()
    // Merge EGM data
    ch_input.map { meta, tsv_list, gbk_list, gbk ->
            return [ meta, tsv_list, gbk ]
        }
        .set { ch_mergemge }
    MERGEEGM( ch_mergemge )

    // Merge GBK files
    MERGEEGM.out.gbk
            .join( ch_input )
            .map { meta, summary_gbk, tsv_list, gbk_list, gbk ->
                return [ meta, summary_gbk, gbk_list ]
            }
            .set { ch_mergegbk }
    MERGEGBKS( ch_mergegbk )

    // Create circos plot
    MERGEGBKS.out.gbk.map
                { meta, name, gbk ->
                return [ meta, name, gbk, "" ]
                }.set { ch_visual }
    VISUALIZE_CIRCULAR(ch_visual)

    emit:
    gbk      = MERGEGBKS.out.gbk     // channel: [meta, markers.gbk]
    tsv      = MERGEEGM.out.tsv      // channel: [meta, summary.tsv]
    png      = VISUALIZE_CIRCULAR.out.png    // channel: [meta, summary.png]
}
