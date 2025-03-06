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
include { COPLA } from '../modules/local/copla/main'
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
        new_name=\$(basename "\$i" | rev | cut -f1 -d'/' | cut -f2- -d'.' | cut -f1 -d'_' | rev)
        echo "Renaming \$i as \${new_name}_${meta.id}.fasta"
        mv "\$i" \${new_name}_${meta.id}.fasta
    done
    """
}

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    // MOBSUITE RECON
    ch_plasmids = MOBSUITE_RECON (ch_samplesheet).plasmids
    ch_versions = ch_versions.mix( MOBSUITE_RECON.out.versions )

    // RENAME PLASMIDS
    // OJO! Ya no filtra por tamaño
    ch_filtered = RENAME_PLASMIDS (ch_plasmids)

    // COPLA
    ch_filtered.flatMap { meta, fileList ->
    // Si fileList no es una lista, lo tokenizamos (dividimos por espacios)
    def files = fileList instanceof List ? fileList : fileList.toString().tokenize(' ')
    files.collect { file -> 
        def plasmid_name = file.getName().toString().replaceFirst(/\.fasta$/, '')
        tuple(meta, plasmid_name, file) }
    } | COPLA

    ch_versions = ch_versions.mix( COPLA.out.versions )

    INTEGRON_FINDER (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix( INTEGRON_FINDER.out.versions )

    // INTEGRON_FINDER.out.view{ it }
    //ch_versions = ch_versions.mix( INTEGRON_FINDER.out.versions )

    // ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // //
    // // Collate and save software versions
    // //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(
    //         storeDir: "${params.outdir}/pipeline_info",
    //         name: 'nf_core_'  +  'pitisfinder_software_'  + 'mqc_'  + 'versions.yml',
    //         sort: true,
    //         newLine: true
    //     ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    emit:
    // multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
