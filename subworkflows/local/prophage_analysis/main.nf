include { PHISPY         } from '../../../modules/nf-core/phispy/main'
include { PROCESS_PHISPY } from '../../../modules/local/processphispy/main'
include { VISUALIZE_LINEAR } from '../../../modules/local/visualize/linear/main'

workflow PROPHAGE_ANALYSIS {

    take:
    ch_gbk          // channel: [ val(meta), [ bam ] ]
    phispy_db       // path: Phispy database (optional)

    main:

    ch_versions = Channel.empty()

    // PHISPY
    ch_phispydb = phispy_db ? Channel.value(phispy_db) : Channel.value([])
    PHISPY(ch_gbk, ch_phispydb)

    ch_versions = ch_versions.mix(PHISPY.out.versions)

    // PROCESS PHISPY
    PROCESS_PHISPY(PHISPY.out.phage_gbk, PHISPY.out.phage_fasta, PHISPY.out.prophage_tsv)

    // CREATE CHANNEL TO LINEAR VISUALIZATION
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

    // VISUALIZATION (LINEAR)
    VISUALIZE_LINEAR ( ch_pplin )

    emit:
    prophages      = PROCESS_PHISPY.out.gbk           // channel: [ val(meta), [ pp_*gbk ] ]
    summary        = PROCESS_PHISPY.out.summary       // channel: [ val(meta), [ prophage_summary.tsv ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

