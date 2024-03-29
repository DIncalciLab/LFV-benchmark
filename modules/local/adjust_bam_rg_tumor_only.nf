process ADJUST_BAM_RG_TUMOR_ONLY {
    tag "Adjust read groups for input BAM: ${meta.sample_name}"
    label 'process_medium'

    container "aldosr/bamsurgeon:1.3-custom"

    input:
    tuple val(meta), val(tumor)
    path picardjar

    output:
    tuple val(meta),    path("*_tumor.bam"),  path("*_tumor.bam.bai"),    emit: tumor_only

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${tumor.tumor_bam} \\
        O=${meta.sample_name}_tumor.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta.sample_name}_tumor \\
        RGLB=${meta.sample_name}_tumor \\
        RGPL=${meta.sample_name}_tumor \\
        RGPU=${meta.sample_name}_tumor \\
        RGSM=${meta.sample_name}_tumor

    samtools index ${meta.sample_name}_tumor.bam
    """

}
