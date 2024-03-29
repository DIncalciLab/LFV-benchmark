process VARSCAN2 {
    tag "Variant calling using VarScan2 on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.4--0':
        'quay.io/biocontainers/varscan:2.4.4--0' }"

    input:
    tuple val(meta), val(normal), val(tumor)
    tuple val(meta), val(mpileup)

    val   fasta
    path  bed

    output:
    //tuple val(meta),
    path("*.vcf"),          optional:true,       emit: vcf_tumor_only
    path("*.snp.vcf.gz"),   optional:true,       emit: vcf_varscan_snp
    path("*.indel.vcf.gz"), optional:true,       emit: vcf_varscan_indel
    //path("*.vcf.tbi")   , emit: vcf_varscan_tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def opt = !params.high_sensitivity ? '' : task.ext.opt
    def opt2 = !params.high_sensitivity ? '' : task.ext.opt2
    def prefix = task.ext.prefix ?: "varscan"
    def mpileup =   !params.tumor_only
                    ? "somatic $mpileup ${prefix} ${opt} --output-vcf --mpileup 1"
                    : "mpileup2cns $mpileup --variants ${opt2} --output-vcf > ${prefix}.vcf"
    def VERSION = '2.4.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    def avail_mem = 3
    if (!task.memory) {
        log.info '[VarScan2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    varscan ${mpileup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """

}