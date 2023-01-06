process VARSCAN2 {
    tag "Variant calling using VarScan2 on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varscan:2.4.4--0':
        'quay.io/biocontainers/varscan:2.4.4--0' }"

    input:
    tuple val(meta), val(tumor_only)
    tuple val(meta), path(normal_bam), path(normal_bai),  path(tumor_bam), path(tumor_bai)
    tuple val(meta), path(mpileup)
    //tuple val(meta), path(bai)

    val   fasta
    path  bed

    output:
    //tuple val(meta),
    path("*.snp.vcf.gz")   , emit: vcf_varscan_snp
    path("*.indel.vcf.gz")   , emit: vcf_varscan_indel
    //path("*.vcf.tbi")   , emit: vcf_varscan_tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "varscan"
    def mpileup = ( normal_bam && tumor_bam )
                    ? "somatic $mpileup ${prefix} --mpileup 1"
                    : "mpileup2cns $mpileup --variants"
    def VERSION = '2.4.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    def avail_mem = 3
    if (!task.memory) {
        log.info '[VarScan2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if ( params.high_sensitivity ) {
    """
    varscan $mpileup \\
        --output-vcf \\
        $args

    gzip -c ${prefix}.snp.vcf > ${prefix}.snp.vcf.gz
    gzip -c ${prefix}.indel.vcf > ${prefix}.indel.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }

    if ( !params.high_sensitivity ){
    """
    varscan $mpileup \\
        --output-vcf

    gzip -c ${prefix}.snp.vcf > ${prefix}.snp.vcf.gz
    gzip -c ${prefix}.indel.vcf > ${prefix}.indel.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }
}