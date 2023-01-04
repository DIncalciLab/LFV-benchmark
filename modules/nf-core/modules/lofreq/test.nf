process TEST {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py36h5b61e8e_8' :
        'quay.io/biocontainers/lofreq:2.1.5--py36h5b61e8e_8' }"

    //container "aldosr/lofreq:2.1.5"

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    val fasta
    path bed

    output:
    path("*.vcf.gz"), emit: test


    script:
    """
    lofreq call -f ${fasta} -o vars.vcf.gz ${normal_bam}

    """}