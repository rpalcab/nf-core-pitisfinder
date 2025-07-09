include { PVOGDOWNLOAD     } from '../../../modules/local/pvogdownload/main'
include { PHISPY           } from '../../../modules/nf-core/phispy/main'
include { PROCESS_PHISPY   } from '../../../modules/local/processphispy/main'
include { VISUALIZE_LINEAR } from '../../../modules/local/visualize/linear/main'
include { GENOMAD_DOWNLOAD } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ANNOTATE } from '../../../modules/local/genomad/annotate/main'
include { GENOMAD_FINDPROVIRUSES } from '../../../modules/local/genomad/findproviruses/main'
include { PROPHAGEMARKERS  } from '../../../modules/local/prophagemarkers/main'
include { PROPHAGEPARSER  } from '../../../modules/local/prophageparser/main'

workflow PROPHAGE_ANALYSIS {

    take:
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    ch_gbk            // channel: [ val(meta), [ bam ] ]
    genomad_db        // path (optional): genomad_db

    main:

    ch_versions = Channel.empty()

    // DOWNLOAD GENOMAD_DB
    ch_genomaddb = Channel.empty()
    if ( !genomad_db ) {
        ch_genomaddb = GENOMAD_DOWNLOAD( ).genomad_db
        ch_versions.mix( GENOMAD_DOWNLOAD.out.versions )
    } else {
        ch_genomaddb = Channel.value(file(genomad_db))
    }

    // GENOMAD
    GENOMAD_ANNOTATE ( ch_fasta, ch_genomaddb )
    ch_fasta.
        join(GENOMAD_ANNOTATE.out.outdir)
        .set { ch_findproviruses }
    GENOMAD_FINDPROVIRUSES ( ch_findproviruses, ch_genomaddb )
    ch_versions.mix( GENOMAD_ANNOTATE.out.versions )

    // PROPHAGE MARKERS
    GENOMAD_FINDPROVIRUSES.out.provirus_genes
            .join ( GENOMAD_FINDPROVIRUSES.out.provirus )
            .join( ch_gbk )
            .set { ch_updategbk }

    PROPHAGEMARKERS ( ch_updategbk )

    // PROPHAGE PARSER
    PROPHAGEMARKERS.out.gbk
            .join( GENOMAD_FINDPROVIRUSES.out.provirus )
            .join( GENOMAD_FINDPROVIRUSES.out.taxonomy )
            .set { ch_prophageparser }

    PROPHAGEPARSER ( ch_prophageparser )

    // // PHISPY
    // PVOGDOWNLOAD()
    // PHISPY(ch_gbk, PVOGDOWNLOAD.out.pvogs_db)

    // ch_versions = ch_versions.mix(PHISPY.out.versions)

    // // PROCESS PHISPY
    // PROCESS_PHISPY(PHISPY.out.phage_gbk, PHISPY.out.phage_fasta, PHISPY.out.prophage_tsv)

    // CREATE CHANNEL TO LINEAR VISUALIZATION
    PROPHAGEPARSER.out.gbk
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
    ch_outvisualize = Channel.value('prophages/summary/')
    VISUALIZE_LINEAR ( ch_pplin, ch_outvisualize )

    emit:
    genomic_gbk    = PROPHAGEMARKERS.out.gbk          // channel: [ val(meta), [ gbk ] ]
    summary        = PROPHAGEPARSER.out.summary      // channel: [ val(meta), [ prophage_summary.tsv ] ]
    gbk            = PROPHAGEPARSER.out.gbk          // channel: [ val(meta), [ prophage_*.gbk ] ]
    png            = VISUALIZE_LINEAR.out.png         // channel: [ val(meta), [ ph_name ] , [ ph_*.gbk ] ]
    versions       = ch_versions                      // channel: [ versions.yml ]
}

