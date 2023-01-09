process ADJUST_BAM_RG {
    tag "Adjust read groups for input BAM: ${meta.sample_name}"
    label 'process_medium'

    container "aldosr/bamsurgeon:1.3-custom"

    input:
    tuple val(meta), val(tumor_only)
    tuple val(meta), val(normal), val(tumor)
    path picardjar

    output:
    tuple val(meta), path("*_tumor.only.bam"),  path("*_tumor.only.bam.bai"),    optional:true,  emit: tumor_only_bam
    tuple val(meta), path("*_normal.bam"), path("*_normal.bam.bai"), path("*_tumor.bam"), path("*_tumor.bam.bai"), optional: true, emit: paired_bam

    when:
    task.ext.when == null || task.ext.when

    script:

    if ( tumor_only ){
    """
        java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${tumor_only.tumor_bam} \
        O=${meta.sample_name}_tumor.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta.sample_name}_tumor \\
        RGLB=${meta.sample_name}_tumor \\
        RGPL=${meta.sample_name}_tumor \\
        RGPU=${meta.sample_name}_tumor \\
        RGSM=${meta.sample_name}_tumor

    samtools index ${meta.sample_name}_tumor.only.bam
    """
    }

    if ( normal && tumor ){
    """
    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${normal.normal_bam} \
        O=${meta.sample_name}_normal.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=${meta.sample_name}_normal \\
        RGLB=${meta.sample_name}_normal \\
        RGPL=${meta.sample_name}_normal \\
        RGPU=${meta.sample_name}_normal \\
        RGSM=${meta.sample_name}_normal

    samtools index ${meta.sample_name}_normal.bam

    java -jar ${picardjar} \\
        AddOrReplaceReadGroups I=${tumor.tumor_bam} \
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
}
