process TEST {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py36h5b61e8e_8' :
        'quay.io/biocontainers/lofreq:2.1.5--py36h5b61e8e_8' }"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    val fasta
    path bed

    output:
    path("*_somatic_final.snvs.vcf.gz"),                                emit: vcf_lofreq_snvs
    path("*_somatic_final_minus-dbsnp.snvs.vcf.gz"),   optional: true,  emit: vcf_lofreq_snvs_minus_dbsnp
    path("*_somatic_final.indels.vcf.gz"),             optional: true,  emit: vcf_lofreq_indels
    path("*_somatic_final_minus-dbsnp.indels.vcf.gz"), optional: true,  emit: vcf_lofreq_indels_minus_dbsnp
    path ("versions.yml"),                                              emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def prefix = task.ext.prefix ?: "lofreq"
    def bam = (normal_bam && tumor_bam) ? "somatic -n ${normal_bam} -t ${tumor_bam} --threads ${task.cpus}" : ''
    """
    lofreq ${bam} -f ${fasta} -o ${suffix}_

    """

    }