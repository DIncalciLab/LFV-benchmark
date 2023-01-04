process TEST {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_low'



    container "aldosr/lofreq:2.1.5"

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