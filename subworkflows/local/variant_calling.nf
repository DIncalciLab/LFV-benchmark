//tag "Variant calling on $bam"

include { GATK4_MUTECT2 }               from '../../modules/nf-core/modules/gatk4/mutect2/main.nf'
include { VARDICTJAVA }                 from '../../modules/nf-core/modules/vardictjava/main.nf'
include { SAMTOOLS_MPILEUP }            from '../../modules/nf-core/modules/samtools/mpileup/main.nf'
include { VARSCAN2 }                    from '../../modules/nf-core/modules/varscan2/main.nf'
include { LOFREQ }                      from '../../modules/nf-core/modules/lofreq/main.nf'
include { STRELKA_SOMATIC }             from '../../modules/nf-core/modules/strelka/somatic/main.nf'
include { FREEBAYES }                   from '../../modules/nf-core/modules/freebayes/main.nf'


workflow VARIANT_CALLING {
    take:
    normal_bam
    tumor_bam          // channel: [ val(meta), bam, bai ]
    fasta
    bed

    main:
    ch_vardict   = Channel.empty()
    ch_mutect    = Channel.empty()
    ch_varscan   = Channel.empty()
    ch_lofreq    = Channel.empty()
    ch_strelka   = Channel.empty()
    ch_freebayes = Channel.empty()

    ch_versions  = Channel.empty()

    VARDICTJAVA(
        normal_bam,
        tumor_bam,
        fasta,
        bed
    )

    GATK4_MUTECT2(
        normal_bam,
        tumor_bam,
        fasta,
        bed
    )

    SAMTOOLS_MPILEUP(
        normal_bam,
        tumor_bam,
        fasta
    )

    VARSCAN2(
        SAMTOOLS_MPILEUP.out.mpileup,
        fasta,
        bed
    )

    LOFREQ(
        normal_bam,
        tumor_bam,
        fasta,
        bed
    )

    STRELKA(
        normal_bam,
        tumor_bam,
        fasta,
        bed
    )

    FREEBAYES(
        normal_bam,
        tumor_bam,
        fasta,
        bed
    )


    ch_vardict   = ch_vardict.mix(VARDICTJAVA.out.vcf_vardict)
    ch_mutect    = ch_mutect.mix(GATK4_MUTECT2.out.vcf_mutect)
    ch_varscan   = ch_varscan.mix(VARSCAN2.out.vcf_varscan)
    ch_lofreq    = ch_vardict.mix(VARDICTJAVA.out.vcf_vardict)
    ch_strelka   = ch_mutect.mix(GATK4_MUTECT2.out.vcf_mutect)
    ch_freebayes = ch_varscan.mix(VARSCAN2.out.vcf_varscan)

    emit:
    vcf_vardict   = ch_vardict       //   channel: [ val(meta), vcf ]
    vcf_mutect    = ch_mutect        //   channel: [ val(meta), vcf ]
    vcf_varscan   = ch_varscan       //   channel: [ val(meta), vcf ]
    vcf_lofreq    = ch_lofreq        //   channel: [ val(meta), vcf ]
    vcf_strelka   = ch_strlka        //   channel: [ val(meta), vcf ]
    vcf_freebayes = ch_freebayes     //   channel: [ val(meta), vcf ]

    //versions  = ch_versions        //   path: versions.yml*/
}
