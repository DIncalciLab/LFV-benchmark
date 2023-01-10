process ADJUST_BAM_RG {
    tag "Adjust read groups for input BAM: ${meta_tumor.sample_name}"
    label 'process_medium'

    container "aldosr/bamsurgeon:1.3-custom"

    input:
    tuple val(meta_tumor), val(tumor)
    tuple val(meta_normal), val(normal)
    path picardjar

    output:
    tuple val(meta_tumor),            path("*_tumor.bam"),  path("*_tumor.bam.bai"),                  emit: tumor_bam
    tuple val(meta_normal),           path("*_normal.bam"), path("*_normal.bam.bai"), optional: true, emit: normal_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${tumor.tumor_bam} \\
        O=${meta_tumor.sample_name}_tumor.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta_tumor.sample_name}_tumor \\
        RGLB=${meta_tumor.sample_name}_tumor \\
        RGPL=${meta_tumor.sample_name}_tumor \\
        RGPU=${meta_tumor.sample_name}_tumor \\
        RGSM=${meta_tumor.sample_name}_tumor

    samtools index ${meta.sample_name}_tumor.bam
    """

    if ( !( normal.normal_bam.isEmpty() ) ) {

    """
    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${normal.normal_bam} \\
        O=${meta_normal.sample_name}_normal.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta_normal.sample_name}_normal \\
        RGLB=${meta_normal.sample_name}_normal \\
        RGPL=${meta_normal.sample_name}_normal \\
        RGPU=${meta_normal.sample_name}_normal \\
        RGSM=${meta_normal.sample_name}_normal

    samtools index ${meta_normal.sample_name}_normal.bam
    """
    }


}
