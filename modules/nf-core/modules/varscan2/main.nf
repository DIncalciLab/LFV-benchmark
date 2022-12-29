process VARSCAN2 {
    tag "Variant calling using VarScan2 on BAMSurgeon spiked-in sample: ${meta.sample}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.4--0':
        'quay.io/biocontainers/varscan:2.4.4--0' }"

    input:
    tuple val(meta), path(mpileup_tumor)
    tuple val(meta), path(mpileup_normal)
    //tuple val(meta), path(bai)

    val   fasta
    path  bed

    output:
    tuple val(meta), path("*.vcf")   , emit: vcf_varscan
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "varscan"
    def tumor_only = ( !mpileup_normal && mpileup_tumor ) ? "mpileup2cns $mpileup_tumor" : ""
    def paired = ( mpileup_normal && mpileup_tumor ) ? "somatic $mpileup_normal $mpileup_tumor" : ""
    def VERSION = '2.4.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    def avail_mem = 3
    if (!task.memory) {
        log.info '[VarScan2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if ( params.high_sensitivity ) {
    """
    varscan $tumor_only \\
        $paired  \\
        $args \\
        --variants > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }

    if ( !params.high_sensitivity ){
    """
    varscan $tumor_only $paired  \\
        --output-vcf \\
        --variants > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }
}