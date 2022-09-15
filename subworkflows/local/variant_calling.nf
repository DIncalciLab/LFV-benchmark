//tag "Variant calling on $bam"

include { GATK4_MUTECT2 }               from '../../modules/nf-core/modules/gatk4/mutect2/main.nf'
include { VARDICTJAVA }                 from '../../modules/nf-core/modules/vardictjava/main.nf'
include { SAMTOOLS_MPILEUP }            from '../../modules/nf-core/modules/samtools/mpileup/main.nf'
include { VARSCAN2 }                    from '../../modules/nf-core/modules/varscan2/main.nf'

workflow VARIANT_CALLING {
    take:
    bam          // channel: [ val(meta), bam, bai ]
    fasta
    bed

    main:
    ch_vardict   = Channel.empty()
    ch_mutect   = Channel.empty()
    ch_varscan   = Channel.empty()

    ch_versions  = Channel.empty()

    VARDICTJAVA(
        bam,
        fasta,
        bed
    )

    GATK4_MUTECT2(
        bam,
        fasta,
        bed
    )

    SAMTOOLS_MPILEUP(
        bam,
        fasta
    )

    VARSCAN2(
        SAMTOOLS_MPILEUP.out.mpileup,
        fasta,
        bed
    )

    ch_vardict = ch_vardict.mix(VARDICTJAVA.out.vcf_vardict)
    ch_mutect = ch_mutect.mix(GATK4_MUTECT2.out.vcf_mutect)
    ch_varscan = ch_varscan.mix(VARSCAN2.out.vcf_varscan)

    ch_vardict.view()
    ch_mutect.view()
    ch_varscan.view()

    emit:
    vcf_vardict   = ch_vardict       //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    vcf_mutect   = ch_mutect       //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    vcf_varscan   = ch_varscan      //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]

    //versions  = ch_versions               //    path: versions.yml*/
}