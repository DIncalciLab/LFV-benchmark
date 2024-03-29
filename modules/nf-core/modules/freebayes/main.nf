process FREEBAYES {
    tag 'Variant calling using FreeBayes on BAMSurgeon spiked-in sample: ${meta.sample_name}'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2' :
        'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2' }"

    input:
    tuple val(meta), val(normal), val(tumor)
    path samples
    path populations
    path cnv
    val fasta
    path target_bed


    output:
    path("*.vcf.gz"), emit: vcf_freebayes
    path("*.vcf.gz.tbi"), emit: vcf_freebayes_tbi
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def opt = task.ext.opt ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    def input            =  !params.tumor_only
                            ? "${tumor.tumor_bam} ${normal.normal_bam}"
                            : "${tumor.tumor_bam}"
    def targets_file     = target_bed     ? "--target ${target_bed}"       : ""
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""

    if ( !params.high_sensitivity ) {
    """
    freebayes \\
        -f $fasta \\
        $targets_file \\
        $samples_file \\
        $populations_file \\
        $cnv_file \\
        $input > ${prefix}.vcf

    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
    } else {
    """
    freebayes \\
        -f ${fasta} \\
        ${targets_file} \\
        ${samples_file} \\
        ${populations_file} \\
        ${cnv_file} \\
        ${opt} \\
        ${input} > ${prefix}.vcf

    bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
    }
}
