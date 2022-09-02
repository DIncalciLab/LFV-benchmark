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
    params.bamsurgeon_path,
    params.picardjar
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


//include { NEAT        } from '../modules/local/neat.nf'
include { BAMSURGEON } from '../modules/local/bamsurgeon.nf'
//include { BENCHMARK  }  from '../modules/local/benchmark.nf'
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
include { GATK4_MUTECT2 }               from '../modules/nf-core/modules/gatk4/mutect2/main'
include { VARDICTJAVA }                 from '../modules/nf-core/modules/vardictjava/main'
include { VARSCAN2 }                    from '../modules/nf-core/modules/varscan2/main'
include { MULTIQC }                     from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow LOWFRAC_VARIANT_BENCHMARK {

    ch_versions = Channel.empty()

    ch_input = Channel
    .fromPath(params.input)
    .splitCsv(header:true, quote:'\"', sep: ",")
    .map { row -> [sample: row.sample, info: row.info]
                }
/*
    ch_rng = Channel
      .from( 1..32767 )
      .randomSample( 1 )
      .view()*/

/*
    NEAT(
        ch_input,
        params.readlen,
        params.coverage,
        params.bed, 
        params.fasta,
        params.neat_path,
        params.fraglen_model,
        params.error_model,
        params.mutation_model,
        params.gc_model
        )
    
    ch_versions = ch_versions.mix(NEAT.out.versions)*/
    
    BAMSURGEON(
        ch_input,
        params.mut_number,
        params.min_fraction,
        params.max_fraction,
        params.maxlen,
        params.fasta,
        params.type,
        params.bed,
        params.picardjar
    )
    
    //ch_versions = ch_versions.mix(BAMSURGEON.out.versions)

    /*
    VARDICTJAVA(
        BAMSURGEON.out.bam,
        //BAMSURGEON.out.bai,
        params.fasta,
        params.bed
    )

    GATK4_MUTECT2(
        BAMSURGEON.out.bam,
        params.bed,
        params.fasta
    )

    VARSCAN2(
        BAMSURGEON.out.bam,
        //BAMSURGEON.out.bai,
        params.bed,
        params.fasta
    )

    /*BENCHMARK(
        VARDICTJAVA.out.vcf, 
        GATK4_MUTECT2.out.vcf,
        VARSCAN2.out.vcf)*/
    
    //ch_versions = ch_versions.mix(BAMSURGEON.out.versions)
    
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    
    workflow_summary    = WorkflowLowFrac.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
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

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
