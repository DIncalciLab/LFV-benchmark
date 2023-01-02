process ADJUST_BAM_RG {
    tag "Adjust read groups for input BAM: ${meta.sample_name}"
    label 'process_low'

    container "aldosr/bamsurgeon:1.3-custom"

    input:
    tuple val(meta), val(normal), val(tumor)
    path picardjar

    output:
    tuple val(meta), val("*_normal.bam"), val("*_normal.bam.bai"), val("*_tumor.bam"), val("*_tumor.bam.bai"), emit: paired_bam

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${normal.normal_bam} \
        O=${meta.sample_name}_normal \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta.sample_name}_normal \\
        RGLB=${meta.sample_name}_normal \\
        RGPL=${meta.sample_name}_normal \\
        RGPU=${meta.sample_name}_normal \\
        RGSM=${meta.sample_name}_normal

    samtools index ${meta.sample_name}_normal.bam

    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${tumor.tumor_bam} \
        O=${meta.sample_name}_tumor \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta.sample_name}_normal \\
        RGLB=${meta.sample_name}_normal \\
        RGPL=${meta.sample_name}_normal \\
        RGPU=${meta.sample_name}_normal \\
        RGSM=${meta.sample_name}_normal

    samtools index ${meta.sample_name}_tumor.bam
    """
}
