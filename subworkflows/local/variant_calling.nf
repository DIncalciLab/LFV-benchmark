/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT CALLING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GATK4_MUTECT2 }               from '../../modules/nf-core/modules/gatk4/mutect2/main.nf'
include { VARDICTJAVA }                 from '../../modules/nf-core/modules/vardictjava/main.nf'
include { SAMTOOLS_MPILEUP }            from '../../modules/nf-core/modules/samtools/mpileup/main.nf'
include { VARSCAN2 }                    from '../../modules/nf-core/modules/varscan2/main.nf'
include { LOFREQ }                      from '../../modules/nf-core/modules/lofreq/main.nf'
include { STRELKA_SOMATIC }             from '../../modules/nf-core/modules/strelka/somatic/main.nf'
include { FREEBAYES }                   from '../../modules/nf-core/modules/freebayes/main.nf'


workflow VARIANT_CALLING {
    take:
    tumor_only
    paired
    germline_resource
    panel_of_normals
    dbsnp_vcf
    manta_candidate_small_indels
    freebayes_samples
    freebayes_population
    freebayes_cnv
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
        tumor_only,
        paired,
        fasta,
        bed
    )

    GATK4_MUTECT2(
        tumor_only,
        paired,
        germline_resource,
        panel_of_normals,
        fasta,
        bed
    )

    SAMTOOLS_MPILEUP(
        tumor_only,
        paired,
        fasta
    )

    VARSCAN2(
        SAMTOOLS_MPILEUP.out.mpileup,
        fasta,
        bed
    )

    LOFREQ (
        tumor_only,
        paired,
        fasta,
        bed,
        dbsnp_vcf
    )

    STRELKA_SOMATIC(
        tumor_only,
        paired,
        manta_candidate_small_indels,
        fasta,
        bed
    )

    FREEBAYES(
        tumor_only,
        paired,
        freebayes_samples,
        freebayes_population,
        freebayes_cnv,
        fasta,
        bed
    )

    ch_output = VARDICTJAVA.out.vcf_vardict
                .mix(   GATK4_MUTECT2.out.vcf_mutect,
                        VARSCAN2.out.vcf_varscan,
                        FREEBAYES.out.vcf_freebayes,
                        LOFREQ.out.vcf_lofreq_snvs,
                        LOFREQ.out.vcf_lofreq_indels,
                        STRELKA_SOMATIC.out.vcf_strelka_snvs,
                        STRELKA_SOMATIC.out.vcf_strelka_indels )

    //ch_vardict   = ch_vardict.mix(VARDICTJAVA.out.vcf_vardict)
    //ch_mutect    = ch_mutect.mix(GATK4_MUTECT2.out.vcf_mutect)
    //ch_varscan   = ch_varscan.mix(VARSCAN2.out.vcf_varscan)
    //ch_lofreq    = ch_vardict.mix(VARDICTJAVA.out.vcf_vardict)
    //ch_strelka   = ch_mutect.mix(GATK4_MUTECT2.out.vcf_mutect)
    //ch_freebayes = ch_varscan.mix(VARSCAN2.out.vcf_varscan)

    emit:
    vcf =  ch_output
    //vcf_vardict   = ch_vardict       //   channel: [ val(meta), vcf ]
    //vcf_mutect    = ch_mutect        //   channel: [ val(meta), vcf ]
    //vcf_varscan   = ch_varscan       //   channel: [ val(meta), vcf ]
    //vcf_lofreq    = ch_lofreq        //   channel: [ val(meta), vcf ]
    //vcf_strelka   = ch_strelka        //   channel: [ val(meta), vcf ]
    //vcf_freebayes = ch_freebayes     //   channel: [ val(meta), vcf ]

    //versions  = ch_versions        //   path: versions.yml*/
}
