process LOFREQ_INDEL {
    tag "Variant calling using LoFreq on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py36h5b61e8e_8' :
        'quay.io/biocontainers/lofreq:2.1.5--py36h5b61e8e_8' }"

    input:
    tuple val(meta), val(normal), val(tumor)
    val fasta
    val bed

    output:
    path("*_somatic_final.snvs.vcf.gz"),               optional: true,  emit: vcf_lofreq_snvs
    path("*_somatic_final_minus-dbsnp.snvs.vcf.gz"),   optional: true,  emit: vcf_lofreq_snvs_minus_dbsnp
    path("*_somatic_final.indels.vcf.gz"),             optional: true,  emit: vcf_lofreq_indels
    path("*_somatic_final_minus-dbsnp.indels.vcf.gz"), optional: true,  emit: vcf_lofreq_indels_minus_dbsnp
    path("*.vcf"),                                     optional:true,   emit: vcf_tumor_only
    path ("versions.yml"),                                              emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def prefix = task.ext.prefix ?: "lofreq"
    def bam = ( !( {assert ${normal.normal_bam} == 'EMPTY'} )  )
                ? "somatic  -n ${prefix}_indel_processed_normal.bam  -t ${prefix}_indel_processed_tumor.bam"
                    "--call-indels -f ${fasta}  -l ${bed}  -o ${prefix}_"
                : "call  --call-indels -f ${fasta}  -l ${bed}  -o ${prefix}.vcf  ${prefix}_indel_processed_tumor.bam"
    def VERSION = '2.1.5'

    def avail_mem = 3
    if (!task.memory) {
        log.info '[LoFreq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (( !( {assert ${normal.normal_bam} == 'EMPTY'} )  )){
    """
    lofreq indelqual \\
        --dindel \\
        -f ${fasta} \\
        -o ${prefix}_indel_processed_normal.bam \\
        ${normal.normal_bam}
    samtools index ${prefix}_indel_processed_normal.bam
    """
    }

    if ( !params.high_sensitivity ) {

    """
    lofreq indelqual \\
        --dindel \\
        -f ${fasta} \\
        -o ${prefix}_indel_processed_tumor.bam \\
        ${tumor.tumor_bam}

    samtools index ${prefix}_indel_processed_tumor.bam

    lofreq ${bam}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq_version: $VERSION
    END_VERSIONS
    """
    }

    }