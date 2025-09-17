include { ABRICATE_RUN as ABRICATE_RUN_NCBI } from '../../../modules/nf-core/abricate/run/main'
include { ABRICATE_RUN as ABRICATE_RUN_VFDB } from '../../../modules/nf-core/abricate/run/main'
include { DEFENSEFINDER_UPDATE              } from '../../../modules/local/defensefinder/update/main'
include { DEFENSEFINDER_RUN                 } from '../../../modules/local/defensefinder/run/main'

workflow RVD_ANNOTATION {

    take:
    ch_fasta    // channel: [ val(meta), [ fasta ] ]
    df_db        // path (optional): df_db

    main:

    ch_versions = Channel.empty()

    // AMR
    ABRICATE_RUN_NCBI ( ch_fasta, [], 'ncbi' )
    ch_versions = ch_versions.mix( ABRICATE_RUN_NCBI.out.versions.first() )

    //VR
    ABRICATE_RUN_VFDB ( ch_fasta, [], 'vfdb' )
    ch_versions = ch_versions.mix( ABRICATE_RUN_VFDB.out.versions.first() )

    //DF
    ch_dfdb = Channel.empty()
    if (!df_db){
        DEFENSEFINDER_UPDATE ()
        ch_dfdb = DEFENSEFINDER_UPDATE.out.db
        ch_versions = ch_versions.mix( DEFENSEFINDER_UPDATE.out.versions )
    } else {
        ch_dfdb = Channel.value(file(df_db))
    }
    DEFENSEFINDER_RUN (
        ch_fasta,
        ch_dfdb
    )
    ch_versions = ch_versions.mix( DEFENSEFINDER_RUN.out.versions.first() )

    emit:
    amr_report      = ABRICATE_RUN_NCBI.out.report      // channel: [ val(meta), [ report ] ]
    vf_report      = ABRICATE_RUN_VFDB.out.report       // channel: [ val(meta), [ report ] ]
    df_report      = DEFENSEFINDER_RUN.out.genes        // channel: [ val(meta), [ report ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

