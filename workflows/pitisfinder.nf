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
include { INTEGRON_FINDER } from '../modules/local/integronfinder/main'
include { INTEGRON_PARSER } from '../modules/local/integronparser/main'
include { IS_BLAST } from '../modules/local/isblast/main'
include { IS_PARSER } from '../modules/local/isparser/main'

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
    ch_samplesheet.map { meta, fasta, ann, res ->
        return [ meta, fasta ]
    }
    .set { ch_fasta }

    //
    // PLASMIDS (ESTO IRÁ A SUBWORKFLOW)
    //

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

    //
    // INTEGRONS (ESTO IRÁ A SUBWORKFLOW)
    //
    // Integron_finder
    INTEGRON_FINDER (ch_fasta)
    ch_versions = ch_versions.mix( INTEGRON_FINDER.out.versions )
    ch_integron_raw = INTEGRON_FINDER.out.integrons
    // Process results
    ch_merged = ch_samplesheet
        .join(ch_integron_raw)
    INTEGRON_PARSER (ch_merged)

    //
    // INSERTION SEQUENCES (ESTO IRÁ A SUBWORKFLOW)
    //
    ch_isdb = Channel.empty()
    if (params.is_db){ 
        ch_isdb = Channel.value(params.is_db)
        // BLASTn search
        IS_BLAST (ch_mobsuite_chr, ch_isdb)
        ch_versions = ch_versions.mix( IS_BLAST.out.versions )
        // Results filtering
        ch_raw_is = IS_BLAST.out.report
        IS_PARSER (ch_raw_is)
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
