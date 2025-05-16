/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { FASTQC                   } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                  } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_pitisfinder_pipeline'
include { MOBSUITE_RECON              } from '../modules/nf-core/mobsuite/recon/main'
include { PLASMIDFINDER               } from '../modules/nf-core/plasmidfinder/main'
include { PLASMID_PARSER               } from '../modules/local/plasmidparser/main'
// include { GENOMAD_DOWNLOAD            } from '../modules/nf-core/genomad/download/main'
// include { GENOMAD_ENDTOEND            } from '../modules/nf-core/genomad/endtoend/main'
include { COPLA_COPLADBDOWNLOAD       } from '../modules/local/copla/copladbdownload/main'
include { COPLA_COPLA                 } from '../modules/local/copla/copla/main'
include { INTEGRONFINDER              } from '../modules/local/integronfinder/main'
include { INTEGRON_PARSER             } from '../modules/local/integronparser/main'
include { IS_BLAST                    } from '../modules/local/isblast/main'
include { IS_PARSER                   } from '../modules/local/isparser/main'
include { ISESCAN                     } from '../modules/local/isescan/main'
include { PHISPY                      } from '../modules/nf-core/phispy/main'
include { PHIGARO_SETUP               } from '../modules/local/phigaro/phigaro_setup'
include { PHIGARO                     } from '../modules/local/phigaro/phigaro'
include { PHASTEST_PHASTESTDBDOWNLOAD } from '../modules/local/phastest/phastestdbdownload/main'
include { PHASTEST_PHASTEST           } from '../modules/local/phastest/phastest/main'
include { MACSYFINDER                 } from '../modules/local/macsyfinder/macsyfinder/main'
include { VERIFYMODEL                 } from '../modules/local/macsyfinder/verifymodel/main'
include { ICEBERG_DB_DOWNLOAD         } from '../modules/local/iceberg/dbdownload/main'
include { ICEBERG_ICESEARCH           } from '../modules/local/iceberg/icesearch/main'
include { ICEBERG_FILTER              } from '../modules/local/iceberg/icefilter/main'
include { ICEFINDER2                  } from '../modules/local/icefinder2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

process FILTER_FAA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(faa), path(gbk), path(chr)

    output:
    tuple val(meta), path("${meta.id}/*_filtered.faa"), emit: chr

    script:
    def prefix = "${meta.id}"
    """
    cut -f2,5 $chr | grep "chromosome" | cut -f2 > chr.txt
    faa_contig_filter.py -f $faa -g $gbk -c chr.txt -o $prefix
    """
}

process MERGE_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(amr), path(gbk)

    output:
    tuple val(meta), path("${meta.id}_merged.gbk"), emit: gbk

    script:
    def prefix = "${meta.id}"
    """
    merge_amr.py -t $amr -g $gbk -o "$prefix"_merged.gbk
    """
}

workflow PITISFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    // Channel solo con sample y fasta
    ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        return [ meta, fasta ]
    }
    .set { ch_fasta }

    ch_samplesheet.map {  meta, fasta, faa, gbk, amr ->
        return [ meta, amr, gbk ]
    }.set {ch_mergeann}
    MERGE_ANNOTATIONS (ch_mergeann)

    ch_samplesheet
        .join(MERGE_ANNOTATIONS.out.gbk)
        .map {  meta, fasta, faa, gbk, amr, merged_gbk ->
            return [ meta, fasta, faa, merged_gbk ]
        }.set {ch_full}

    //
    // PLASMIDS (ESTO IRÁ A SUBWORKFLOW)
    //
    if ( !params.skip_plasmids ) {
        //
        // MOBSUITE RECON
        //
        MOBSUITE_RECON (ch_fasta)
        ch_mobsuite_pl = MOBSUITE_RECON.out.plasmids
        ch_mobsuite_chr = MOBSUITE_RECON.out.chromosome
        ch_versions = ch_versions.mix( MOBSUITE_RECON.out.versions )
        //
        // PLASMIDFINDER
        //
        // PLASMIDFINDER (ch_fasta)
        // ch_versions = ch_versions.mix( PLASMIDFINDER.out.versions )
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

        // PREPARE PLASMID_PARSER CHANNEL
        COPLA_COPLA.out.query
            .join(COPLA_COPLA.out.ptu, by: [0, 1])
            .set { ch_coplajoint }

        ch_full.map { meta, fasta, faa, gbk ->
            return [ meta, gbk ]
            }
            .join(MOBSUITE_RECON.out.mobtyper_results)
            .join(MOBSUITE_RECON.out.contig_report)
            .set { ch_mobsample }

        ch_mobsample
            .cross(ch_coplajoint)
            .map { mob, copla ->
                def meta = mob[0]
                def gbk = mob[1]
                def mob_typer = mob[2]
                def contig_report = mob[3]
                def plasmid_name = copla[1]
                def qry = copla[2]
                def ptu = copla[3]
                return [ meta, plasmid_name, qry, ptu, gbk, mob_typer, contig_report ]
            }
            .set { ch_plasmidparser }

        PLASMID_PARSER (ch_plasmidparser)
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
        ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
            return [ meta, fasta, amr ]
            }.join(ch_integron_raw)
            .set { ch_merged }
        // INTEGRON_PARSER (ch_merged)
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
        // ch_isdb = Channel.empty()
        // if (params.is_db){
        //     ch_isdb = Channel.value(params.is_db)
        //     // BLASTn search
        //     IS_BLAST (ch_is_input, ch_isdb)
        //     ch_versions = ch_versions.mix( IS_BLAST.out.versions )
        //     // Results filtering
        //     ch_raw_is = IS_BLAST.out.report
        //     IS_PARSER (ch_raw_is)
        // }
        ISESCAN (ch_is_input)
        ch_versions = ch_versions.mix( ISESCAN.out.versions )
    }

    if ( !params.skip_prophages ) {
        //
        // PHIPSY
        //
        // ch_phispydb = Channel.value([])
        // if (params.phispy_db){
        //     ch_phispydb = Channel.value(params.phispy_db)
        // }
        // ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        //     return [ meta, gbk ]
        //     }.set { ch_phispy }
        // PHISPY (ch_phispy, ch_phispydb)
        // ch_versions = ch_versions.mix( PHISPY.out.versions )

        //
        // PHIGARO
        //
        ch_phigaro_setup = params.phigaro_db ? Channel.fromPath(params.phigaro_db).map{ it -> [it] } : Channel.of([])
        PHIGARO_SETUP(ch_phigaro_setup)
        // Create value channels for db and config
        ch_db = PHIGARO_SETUP.out.pvog.first()
        ch_config = PHIGARO_SETUP.out.config.first()

        PHIGARO(
            ch_fasta,
            ch_db,
            ch_config
        )
        ch_versions = ch_versions.mix( PHIGARO.out.versions )

        //
        // PHASTEST
        //
        // ch_phastestdb = Channel.empty()
        // if (!params.phastest_db){
        //     PHASTEST_PHASTESTDBDOWNLOAD ()
        //     ch_phastestdb = PHASTEST_PHASTESTDBDOWNLOAD.out.db
        // } else {
        //     ch_phastestdb = Channel.value(params.phastest_db)
        // }
        // PHASTEST_PHASTEST ( ch_fasta, ch_phastestdb)
        // ch_versions = ch_versions.mix( PHASTEST_PHASTEST.out.versions )

        //
        // Genomad
        //

        // if ( params.genomad_db ) {
        //     ch_db_for_genomad = Channel.fromPath(params.genomad_db)
        // } else {
        //     ch_db_for_genomad = GENOMAD_DOWNLOAD( ).genomad_db
        //     ch_versions.mix( GENOMAD_DOWNLOAD.out.versions )
        // }
        // ch_identified_viruses = GENOMAD_ENDTOEND ( ch_fasta, ch_db_for_genomad ).virus_fasta
        // ch_versions.mix( GENOMAD_ENDTOEND.out.versions )
    }

    if ( !params.skip_ices ) {
        //
        // ICEs (ESTO IRÁ A SUBWORKFLOW)
        //
        // MACSYFINDER
        ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
            return [ meta, faa, gbk ]
        }.join(MOBSUITE_RECON.out.contig_report)
        .set { ch_msy_preproc }
        FILTER_FAA (ch_msy_preproc)
        ch_msymodel = Channel.value('CONJScan/Plasmids')
        VERIFYMODEL (ch_msymodel)
        ch_model = VERIFYMODEL.out.model
        MACSYFINDER (FILTER_FAA.out.chr, ch_model)
        ch_versions = ch_versions.mix( MACSYFINDER.out.versions )
        //
        // ICEBERG
        //
        // Create a channel with the ICEberg database path or 'null' if not provided
        // if ( !params.iceberg_db){
        //     // Run the ICEBERG_DB_DOWNLOAD process
        //     ICEBERG_DB_DOWNLOAD()
        //     ch_iceberg_db = ICEBERG_DB_DOWNLOAD.out.db
        // }
        // else {
        //     ch_iceberg_db = Channel.value(params.iceberg_db)
        // }
        // ICEBERG_ICESEARCH(ch_mobsuite_chr, ch_iceberg_db)
        // ch_versions = ch_versions.mix( ICEBERG_ICESEARCH.out.versions )
        // ICEBERG_FILTER(ICEBERG_ICESEARCH.out.tsv)
        //
        // ICEFINDER2
        //
        // ch_samplesheet.map { meta, fasta, faa, gbk, amr ->
        //     return [ meta, gbk ]
        // }.set { ch_icefinder }
        // ICEFINDER2 ( ch_icefinder )
        // ch_versions = ch_versions.mix( ICEFINDER2.out.versions )
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
