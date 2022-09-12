//
// Variant calling
//

include { GATK4_MUTECT2 }               from '../modules/nf-core/modules/gatk4/mutect2/main.nf'
include { VARDICTJAVA }                 from '../modules/nf-core/modules/vardictjava/main.nf'
include { SAMTOOLS_MPILEUP }            from '../modules/nf-core/modules/samtools/mpileup/main.nf'
include { VARSCAN2 }                    from '../modules/nf-core/modules/varscan2/main.nf'

workflow SPIKEIN_SOMATIC {
    take:
    bam          // channel: [ val(meta), bam ]
    bai
    fasta
    bed

    main:
    ch_vardict   = Channel.empty()
    ch_mutect   = Channel.empty()
    ch_varscan   = Channel.empty()

    ch_versions  = Channel.empty()

    VARDICTJAVA(
        bam,
        bai,
        fasta,
        bed
    )

    GATK4_MUTECT2(
        bam,
        bai,
        fasta,
        bed
    )

    SAMTOOLS_MPILEUP(
        bam,
        fasta
    )

    VARSCAN2(
        SAMTOOLS.MPILEUP.out.bam,
        fasta,
        bed
    )

    ch_vardict = ch_vardict.mix(VARDICTJAVA.out.vcf)
    ch_mutect = ch_mutect.mix(GATK4_MUTECT2.out.vcf)
    ch_varscan = ch_varscan.mix(VARSCAN2.out.vcf)

    emit:
    vcf_vardict   = ch_vardict       //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    vcf_mutect   = ch_mutect       //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    vcf_varscan   = ch_varscan      //   channel: [ val(meta), vcf.gz, vcf.gz.tbi ]

    //versions  = ch_versions               //    path: versions.yml
}