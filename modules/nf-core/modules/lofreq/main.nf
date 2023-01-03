process LOFREQ {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39hf2bf078_8':
        'quay.io/biocontainers/lofreq:2.1.5--py39hf2bf078_8' }"

    input:
    tuple val(meta), path(tumor_only)
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)

    val   fasta

    path  bed
    path  dbsnp_vcf

    output:
    //tuple val(meta),
    path("*_somatic_final.snvs.vcf.gz")   , emit: vcf_lofreq_snvs
    path("*_somatic_final_minus-dbsnp.snvs.vcf.gz")   , emit: vcf_lofreq_snvs_minus_dbsnp, optional: true
    path("*_somatic_final.indels.vcf.gz"), emit: vcf_lofreq_indels, optional: true
    path("*_somatic_final_minus-dbsnp.indels.vcf.gz"), emit: vcf_lofreq_indels_minus_dbsnp, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "lofreq"
    def bam    = ( normal_bam && tumor_bam )
                    ? "somatic -n ${normal_bam} -t ${tumor_bam} --threads $task.cpus"
                    : "call-parallel ${tumor_only.tumor_bam} --pp-threads $task.cpus"
    def dbsnp =  dbsnp_vcf ? "--d $dbsnp_vcf" : ""
    def VERSION = '2.1.5'

    def avail_mem = 3
    if (!task.memory) {
        log.info '[LoFreq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (params.high_sensitivity) {
    """
    lofreq $bam \\
        -f $fasta \\
        -o ${prefix}_

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }

    if ( !params.high_sensitivity ){
    """
    lofreq $bam \\
        -f $fasta \\
        -o ${prefix}_

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan2_version: $VERSION
    END_VERSIONS
    """
    }
}