process LOFREQ {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39hf2bf078_8':
        'quay.io/biocontainers/lofreq:2.1.5--py39hf2bf078_8' }"

    input:
    tuple val(meta), path(normal_bam)
    tuple val(meta), path(tumor_bam)

    val   fasta
    path  bed

    output:
    tuple val(meta), path("*.vcf")   , emit: vcf_lofreq
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "lofreq"
    def VERSION = '2.1.5'

    def avail_mem = 3
    if (!task.memory) {
        log.info '[LoFreq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (params.mode == 'high-sensitivity') {
    """
    varscan mpileup2cns $mpileup  \\
        $args \\
        --variants > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    } else {
    """
    varscan mpileup2cns $mpileup  \\
        --output-vcf \\
        --variants > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }
}