include { MACSYFINDER       } from '../../../modules/local/macsyfinder/macsyfinder/main'
include { VERIFYMODEL       } from '../../../modules/local/macsyfinder/verifymodel/main'
include { ICEFINDER2        } from '../../../modules/local/icefinder2/main'
include { FILTER_FAA        } from '../../../modules/local/filterfaa/main'

workflow ICE_ANALYSIS {

    take:
    ch_faa          // channel: [ val(meta), [ faa ] ]
    ch_gbk          // channel: [ val(meta), [ bam ] ]
    ch_contig_report // channel: [ val(meta), [ contig_report ] ]

    main:

    ch_versions = Channel.empty()

    ch_faa
        .join(ch_gbk)
        .join(ch_contig_report)
        .set { ch_msy_preproc }

    FILTER_FAA(ch_msy_preproc)

    // MACSYFINDER
    ch_msymodel = Channel.value('CONJScan/Plasmids')
    VERIFYMODEL(ch_msymodel)
    MACSYFINDER(FILTER_FAA.out.chr, VERIFYMODEL.out.model)
    ch_versions = ch_versions.mix(MACSYFINDER.out.versions)

    // ICEFINDER2
    // ICEFINDER2 ( ch_icefinder )
    // ch_versions = ch_versions.mix( ICEFINDER2.out.versions )

    emit:
    // results         = MACSYFINDER.out

    versions = ch_versions                     // channel: [ versions.yml ]
}

