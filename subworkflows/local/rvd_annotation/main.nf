include { ABRICATE_RUN as ABRICATE_RUN_NCBI } from '../../../modules/nf-core/abricate/run/main'
include { ABRICATE_RUN as ABRICATE_RUN_VFDB } from '../../../modules/nf-core/abricate/run/main'
include { AMRFINDERPLUS_UPDATE              } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                 } from '../../../modules/nf-core/amrfinderplus/run/main'
include { DEFENSEFINDER_UPDATE              } from '../../../modules/local/defensefinder/update/main'
include { DEFENSEFINDER_RUN                 } from '../../../modules/local/defensefinder/run/main'
include { FORMATAMRFINDER                   } from '../../../modules/local/formatamrfinder/main'

workflow RVD_ANNOTATION {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    df_db        // path (optional): df_db
    amr_db       // path (optional): amr_db

    main:

    ch_versions = Channel.empty()

    // AMR
    ch_amr = Channel.empty()
    if ( params.amr_annotator == 'amrfinder' ) {
        if ( params.amr_db ) {
            ch_amrfinder_db = Channel
                .fromPath(params.amr_db, checkIfExists: true)
                .first()
        }
        else if ( !params.amr_db ) {
            AMRFINDERPLUS_UPDATE()
            ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
            ch_amrfinder_db = AMRFINDERPLUS_UPDATE.out.db
        }
        AMRFINDERPLUS_RUN(ch_fasta, ch_amrfinder_db)

        FORMATAMRFINDER( AMRFINDERPLUS_RUN.out.report )
        ch_amr = ch_amr.mix( FORMATAMRFINDER.out.report )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)
    }
    else if ( params.amr_annotator == 'abricate' ) {
        ABRICATE_RUN_NCBI ( ch_fasta, [], 'ncbi' )
        ch_amr = ch_amr.mix( ABRICATE_RUN_NCBI.out.report )
        ch_versions = ch_versions.mix( ABRICATE_RUN_NCBI.out.versions.first() )
    }

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
    amr_report     = ch_amr                             // channel: [ val(meta), [ report ] ]
    vf_report      = ABRICATE_RUN_VFDB.out.report       // channel: [ val(meta), [ report ] ]
    df_report      = DEFENSEFINDER_RUN.out.genes        // channel: [ val(meta), [ report ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

