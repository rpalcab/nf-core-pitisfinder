include { INTEGRONFINDER   } from '../../../modules/local/integronfinder/main'
include { INTEGRONMARKERS  } from '../../../modules/local/integronmarkers/main'
include { INTEGRON_PARSER  } from '../../../modules/local/integronparser/main'
include { VISUALIZE_LINEAR } from '../../../modules/local/visualize/linear/main'

workflow INTEGRON_ANALYSIS {

    take:
    ch_fasta        // channel: [ val(meta), fasta ]
    ch_gbk          // channel: [ val(meta), gbk ]

    main:

    ch_versions = Channel.empty()

    // INTEGRONFINDER
    INTEGRONFINDER ( ch_fasta )
    ch_versions = ch_versions.mix(INTEGRONFINDER.out.versions.first())

    // PREPARE INTEGRONARKERS CHANNEL
    INTEGRONFINDER.out.integrons
        .join( ch_gbk )
        .set { ch_updategbk }

    // INTEGRON MARKERS
    INTEGRONMARKERS( ch_updategbk )

    INTEGRONMARKERS.out.gbk
        .join( INTEGRONFINDER.out.integrons )
        .set { ch_integronparser }

    // INTEGRON_PARSER
    INTEGRON_PARSER ( ch_integronparser )

    // PREPARE VISUALIZATION CHANNEL
    INTEGRON_PARSER.out.gbk
        .flatMap { meta, gbk_files ->
        if (gbk_files instanceof List) {
            // If gbk_files is a list, process each file
            return gbk_files.collect { gbk_file ->
                def int_meta = [id: gbk_file.name.tokenize('.')[0]]
                [meta, int_meta, gbk_file]
            }
        } else {
            // If gbk_files is a single file, process it directly
            def int_meta = [id: gbk_files.name.tokenize('.')[0]]
            return [[meta, int_meta, gbk_files]]
        }
    }.set { ch_vislin }

    // VISUALIZATION (LINEAR)
    VISUALIZE_LINEAR ( ch_vislin )

    emit:
    integrons      = INTEGRONFINDER.out.integrons     // channel: [ val(meta), [ *.integrons ] ]
    summary        = INTEGRON_PARSER.out.summary      // channel: [ val(meta), [ integron_summary.tsv ] ]
    gbk            = INTEGRON_PARSER.out.gbk          // channel: [ val(meta), [ int_*.gbk ] ]
    genomic_gbk    = INTEGRONMARKERS.out.gbk          // channel: [ val(meta), [ gbk ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

