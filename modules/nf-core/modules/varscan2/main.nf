process VARSCAN2 {
    //tag "$meta.id"
    label 'process_medium'

    // ADJUST THIS SECTION
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::varscan=2.4.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0':
        'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0' }"

    input:
    //tuple val(meta), path(bam), path(bai), path(bed)
    //tuple path(fasta), path(fasta_fai)
    path bam
    path bed
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf_varscan
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "varscan"
    def VERSION = '2.4.4' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    
    """
    samtools mpileup \\
        $args \\
        -f $fasta \\
        "${suffix}.bam" |
                varscan mpileup2cns  \\
                    $args \\
                    --variants > "${prefix}.vcf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
}