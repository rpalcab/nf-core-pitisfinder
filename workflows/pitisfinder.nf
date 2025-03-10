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
    ch_mobsuite = MOBSUITE_RECON (ch_fasta).plasmids
    ch_versions = ch_versions.mix( MOBSUITE_RECON.out.versions )

    // RENAME PLASMIDS
    // OJO! Ya no filtra por tamaño
    ch_rename = RENAME_PLASMIDS (ch_mobsuite)

    
    ch_rename.flatMap { meta, fileList ->
        // Si fileList no es una lista, lo tokenizamos (dividimos por espacios)
        def files = fileList instanceof List ? fileList : fileList.toString().tokenize(' ')
        files.collect { file -> 
            def plasmid_name = file.getName().toString().replaceFirst(/\.fasta$/, '')
            tuple(meta, plasmid_name, file) }
    }.set { ch_plasmids }

    // COPLA
    if (!params.copla_db){ 
        COPLA_COPLADBDOWNLOAD ()
        copla_db_path = COPLA_COPLADBDOWNLOAD.out.db
    } else {
        copla_db_path = params.copla_db
    }

    COPLA_COPLA ( ch_plasmids, copla_db_path )

    ch_versions = ch_versions.mix( COPLA_COPLA.out.versions )

    //
    // INTEGRONS (ESTO IRÁ A SUBWORKFLOW)
    //
    INTEGRON_FINDER (ch_fasta)
    ch_versions = ch_versions.mix( INTEGRON_FINDER.out.versions )

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
