/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowLowFrac.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.fasta,
    params.bed,
    params.neat_path,
    params.readlen,
    params.coverage,
    params.error_model,
    params.mutation_model,
    params.gc_model,
    params.fraglen_model,
    params.picardjar,
    params.skip_normal_generation,
    params.skip_tumor_generation,
    params.skip_variant_calling
]
/*
for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
    }
}*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { NEAT        }               from '../modules/local/neat.nf'
include { BAMSURGEON }                from '../modules/local/bamsurgeon.nf'
include { ADJUST_BAM_RG_PAIRED }      from '../modules/local/adjust_bam_rg_paired.nf'
include { ADJUST_BAM_RG_TUMOR_ONLY }  from '../modules/local/adjust_bam_rg_tumor_only.nf'
include { VARIANT_CALLING }           from '../subworkflows/local/variant_calling.nf'
include { GENERATE_PLOTS  }           from '../modules/local/generate_plots.nf'
//include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC }                     from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE CHANNELS BASED ON PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
input_all          = !(params.skip_normal_generation)
                     ? Channel
                     .fromPath(params.input_all)
                     .splitCsv(header:true, quote:'\"', sep: ",")
                     .map { row -> [sample_name: row.sample, info: row.info] }
                     : Channel.value([])

input_normal       = ( params.skip_normal_generation )
                     ? Channel
                     .fromFilePairs(params.input_normal + "/*.{bam,bai}", flat:true )
                     { sample_name -> sample_name.name.replaceAll(/.normal|.bam|.bai$/,'') }
                     .map { sample_name, bam, bed -> [[sample_name: sample_name], [normal_bam: bam, normal_bai: bed ]]}
                     : Channel.empty()

input_tumor        = ( params.skip_normal_generation && params.skip_tumor_generation )
                     ? Channel
                     .fromFilePairs(params.input_tumor + "/*.{bam,bai}", flat:true )
                     { sample_name -> sample_name.name.replaceAll(/.tumor|.bam|.bai$/,'') }
                     .map { sample_name, bam, bed -> [[sample_name: sample_name], [tumor_bam: bam, tumor_bai: bed ]]}
                     : Channel.empty()

germline_resource  = params.germline_resource
                     ? Channel
                     .fromPath(params.germline_resource)
                     .collect()
                     : Channel.value([])

panel_of_normals   = params.panel_of_normals
                     ? Channel
                     .fromPath(params.panel_of_normals)
                     .collect()
                     : Channel.value([])

dbsnp_vcf          = params.dbsnp_vcf          ?: Channel.value([])

manta_candidate_small_indels  = params.manta_candidate_small_indels
                     ? Channel
                     .fromPath(params.manta_candidate_small_indels)
                     .collect()
                     : Channel.value([])

freebayes_samples = params.freebayes_samples
                     ? Channel
                     .fromPath(params.freebayes_samples)
                     .collect()
                     : Channel.value([])

freebayes_population = params.freebayes_population
                     ? Channel
                     .fromPath(params.freebayes_population)
                     .collect()
                     : Channel.value([])

freebayes_cnv = params.freebayes_cnv
                     ? Channel
                     .fromPath(params.freebayes_cnv)
                     .collect()
                     : Channel.value([])


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow LOWFRAC_VARIANT_BENCHMARK {

    ch_versions = Channel.empty()

    if ( !(params.skip_normal_generation) ){

        NEAT(
            input_all,
            params.bed,
            params.fasta
        )
     }
    //ch_versions = ch_versions.mix(NEAT.out.versions)

    if ( !(params.skip_normal_generation) && !(params.skip_tumor_generation) ){

        input_normal = NEAT
                        .out
                        .bam
                        .map { meta, bam, bai ->
                               [
                                    [sample_name: meta.sample_name], [normal_bam: bam, normal_bai: bai ]]}

        BAMSURGEON(
            input_normal,
            params.mut_number,
            params.min_fraction,
            params.max_fraction,
            params.maxlen,
            params.fasta,
            params.bed,
            params.picardjar
        )

        input_tumor = BAMSURGEON
                        .out
                        .bam
                        .map { meta, bam, bai ->
                               [
                                    [sample_name: meta.sample_name], [tumor_bam: bam, tumor_bai: bai ]
                                    ]
                                    }


    }
    //ch_versions = ch_versions.mix(BAMSURGEON.out.versions)

    if ( params.skip_normal_generation && !(params.skip_tumor_generation )){

        BAMSURGEON(
            input_normal,
            params.mut_number,
            params.min_fraction,
            params.max_fraction,
            params.maxlen,
            params.fasta,
            params.bed,
            params.picardjar
        )

        input_tumor = BAMSURGEON
                        .out
                        .bam
                        .map { meta, bam, bai ->
                               [
                                    [sample_name: meta.sample_name], [tumor_bam: bam, tumor_bai: bai ]
                                    ]
                                    }
    }

    if ( !params.skip_variant_calling ){

        if ( params.tumor_only ){

            ADJUST_BAM_RG_TUMOR_ONLY(
                input_tumor,
                params.picardjar
            )
            tumor_adjusted = ADJUST_BAM_RG_TUMOR_ONLY
                             .out
                             .tumor_only
                             .map { it ->
                                [
                                  [sample_name: it[0].sample_name],
                                  [normal_bam: 'EMPTY', normal_bai: 'EMPTY'],
                                  [tumor_bam: it[1],      tumor_bai: it[1]]
                                ]
                                    }
             input_calling = tumor_adjusted
        }
        if (!params.tumor_only) {

            input_paired = (input_normal.join(input_tumor))

            ADJUST_BAM_RG_PAIRED(
                input_paired,
                params.picardjar
            )

            input_calling  = ADJUST_BAM_RG_PAIRED
                              .out
                              .paired
                              .map{ it ->
                                    [
                                        [sample_name: it[0].sample_name],
                                        [normal_bam: it[1], normal_bai: it[2]],
                                        [tumor_bam: it[3], tumor_bai: it[4]]
                                    ]
                                   }
        }

        VARIANT_CALLING(
            input_calling,
            germline_resource,
            panel_of_normals,
            manta_candidate_small_indels,
            freebayes_samples,
            freebayes_population,
            freebayes_cnv,
            dbsnp_vcf,
            params.fasta,
            params.bed
        )
    }
/*

        GENERATE_PLOTS(
            neat_ch.vcf,
            bamsurgeon_ch.vcf,

            vardict_ch.vcf,
            mutect_ch.vcf,
            varscan_ch.vcf
        )
    }

/*if (params.skip_tumor_generation && !(params.skip_variant_calling)){

        VARIANT_CALLING(
            ch_tumor_bam,
            params.fasta,
            params.bed
        )

        GENERATE_PLOTS(
            neat_ch.vcf,
            bamsurgeon_ch.vcf,
            vardict_ch.vcf,
            mutect_ch.vcf,
            varscan_ch.vcf
        )
    }*/

    /*

    if (params.skip_normal_generation && params.skip_tumor_generation && params.skip_variant_calling){
        log.error "You need to specify an option for the pipeline. See the README for help."
        exit 1
    }
    //ch_versions = ch_versions.mix(BAMSURGEON.out.versions)

    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowLowFrac.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    /*
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)*/


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
