process LOFREQ_SNV {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py36h5b61e8e_8' :
        'quay.io/biocontainers/lofreq:2.1.5--py36h5b61e8e_8' }"

    input:
    tuple val(meta), val(normal), val(tumor)
    val fasta
    val bed

    output:
    path("*.vcf"),                                     optional:true,   emit: vcf_tumor_only
    path("*_somatic_final.snvs.vcf.gz"),               optional:true,   emit: vcf_lofreq_snvs
    path("*_somatic_final_minus-dbsnp.snvs.vcf.gz"),   optional: true,  emit: vcf_lofreq_snvs_minus_dbsnp
    path ("versions.yml"),                                              emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def prefix = task.ext.prefix ?: "lofreq"
    def bam = ( !( params.tumor_only )  )
                ? "somatic  -n ${normal.normal_bam}  -t ${tumor.tumor_bam}  -f ${fasta}  -l ${bed}  -o ${prefix}_"
                : "call  -f ${fasta}  -l ${bed}  -o ${prefix}.vcf  ${tumor.tumor_bam}"
    def VERSION = '2.1.5'

    def avail_mem = 3
    if (!task.memory) {
        log.info '[LoFreq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if ( !params.high_sensitivity ) {
    """
    lofreq ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq_version: $VERSION
    END_VERSIONS
    """
    }

    }