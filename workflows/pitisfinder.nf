/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { FASTQC                 } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'
include { MOBSUITE_RECON  } from '../modules/nf-core/mobsuite/recon/main'
include { COPLA_COPLADBDOWNLOAD } from '../modules/local/copla/copladbdownload/main'
include { COPLA_COPLA } from '../modules/local/copla/copla/main'
include { INTEGRONFINDER } from '../modules/local/integronfinder/main'
include { INTEGRON_PARSER } from '../modules/local/integronparser/main'
include { IS_BLAST } from '../modules/local/isblast/main'
include { IS_PARSER } from '../modules/local/isparser/main'
include { ISESCAN } from '../modules/local/isescan/main'
include { PHASTEST_PHASTESTDBDOWNLOAD } from '../modules/local/phastest/phastestdbdownload/main'
include { PHASTEST_PHASTEST } from '../modules/local/phastest/phastest/main'
include { MACSYFINDER } from '../modules/local/macsyfinder/macsyfinder/main'
include { VERIFYMODEL } from '../modules/local/macsyfinder/verifymodel/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ECHO {
    input:
    tuple val(meta), path(fasta)

    output:
    stdout

    script:
    """
    echo $meta $fasta
    """
}

process RENAME_PLASMIDS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(plasmid_files)

    output:
    tuple val(meta), path("*.fasta"), emit: filtered_plasmids, optional: true

    script:
    """
    for i in ${plasmid_files}; do
        new_name=\$(basename "\$i" | sed 's/plasmid_//' | sed 's/.fasta//')
        mv "\$i" \${new_name}_${meta.id}.fasta
    done
    """
}

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    // Channel solo con sample y fasta
    ch_samplesheet.map { meta, fasta, gene_ann, prot_ann, res ->
        return [ meta, fasta ]
    }
    .set { ch_fasta }

    //
    // PLASMIDS (ESTO IRÁ A SUBWORKFLOW)
    //
    if ( !params.skip_plasmids ) {
        // MOBSUITE RECON
        MOBSUITE_RECON (ch_fasta)
        ch_mobsuite_pl = MOBSUITE_RECON.out.plasmids
        ch_mobsuite_chr = MOBSUITE_RECON.out.chromosome
        ch_versions = ch_versions.mix( MOBSUITE_RECON.out.versions )

        // RENAME PLASMIDS
        // OJO! Ya no filtra por tamaño
        ch_rename = RENAME_PLASMIDS (ch_mobsuite_pl)
        ch_rename
            .flatMap { meta, plasmid_files ->
                def file_list = plasmid_files instanceof List ? plasmid_files : [plasmid_files]
                file_list.collect { file ->
                    def plasmid_name = file.baseName
                    return tuple(meta, plasmid_name, file)
                }
            }
            .set{ ch_plasmids }
        // COPLA
        ch_copladb = Channel.empty()
        if (!params.copla_db){
            COPLA_COPLADBDOWNLOAD ()
            ch_copladb = COPLA_COPLADBDOWNLOAD.out.db
        } else {
            ch_copladb = Channel.value(params.copla_db)
        }
        COPLA_COPLA ( ch_plasmids, ch_copladb)
        ch_versions = ch_versions.mix( COPLA_COPLA.out.versions )
    }

    if ( !params.skip_integrons ) {
        //
        // INTEGRONS (ESTO IRÁ A SUBWORKFLOW)
        //
        // INTEGRONFINDER
        INTEGRONFINDER (ch_fasta)
        ch_versions = ch_versions.mix( INTEGRONFINDER.out.versions )
        ch_integron_raw = INTEGRONFINDER.out.integrons
        // Process results
        ch_merged = ch_samplesheet
            .join(ch_integron_raw)
        INTEGRON_PARSER (ch_merged)
    }

    if ( !params.skip_is ) {
        //
        // INSERTION SEQUENCES (ESTO IRÁ A SUBWORKFLOW)
        //
        ch_is_input = Channel.empty()
        if ( !params.skip_plasmids ) {
            ch_is_input = ch_mobsuite_chr
        } else {
            ch_is_input = ch_fasta
        }
        ch_isdb = Channel.empty()
        if (params.is_db){
            ch_isdb = Channel.value(params.is_db)
            // BLASTn search
            IS_BLAST (ch_is_input, ch_isdb)
            ch_versions = ch_versions.mix( IS_BLAST.out.versions )
            // Results filtering
            ch_raw_is = IS_BLAST.out.report
            IS_PARSER (ch_raw_is)
        }
        ISESCAN (ch_is_input)
        ch_versions = ch_versions.mix( ISESCAN.out.versions )
    }

    if ( !params.skip_prophages ) {
        //
        // PHASTEST
        //
        ch_phastestdb = Channel.empty()
        if (!params.phastest_db){
            PHASTEST_PHASTESTDBDOWNLOAD ()
            ch_phastestdb = PHASTEST_PHASTESTDBDOWNLOAD.out.db
        } else {
            ch_phastestdb = Channel.value(params.phastest_db)
        }
        PHASTEST_PHASTEST ( ch_fasta, ch_phastestdb)
        // ch_versions = ch_versions.mix( PHASTEST_PHASTEST.out.versions )
    }

    if ( !params.skip_ices ) {
        //
        // ICEs (ESTO IRÁ A SUBWORKFLOW)
        //
        // MACSYFINDER
        ch_samplesheet.map { meta, fasta, gene_ann, prot_ann, res ->
            return [ meta, prot_ann ]
        }
        .set { ch_macsyfinder }
        ch_msymodel = Channel.value(['CONJScan','CONJScan/Plasmids'])
        VERIFYMODEL (ch_msymodel)
        ch_model = VERIFYMODEL.out.model
        MACSYFINDER (ch_macsyfinder, ch_model)
        ch_versions = ch_versions.mix( MACSYFINDER.out.versions )
    }

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
