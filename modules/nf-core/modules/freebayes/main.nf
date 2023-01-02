process FREEBAYES {
    tag 'Variant calling using FreeBayes on BAMSurgeon spiked-in sample: ${meta.sample_name}'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2' :
        'quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2' }"

    input:
    tuple val(meta), val(tumor_only)
    tuple val(meta), path(normal_bam), path(normal_bai), path(tumor_bam), path(tumor_bai)
    path samples
    path populations
    path cnv
    path fasta
    path target_bed
    //path fasta_fai


    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.isample_name}"
    def input            = (tumor && normal)       ? "${tumor_bam} ${normal_bam}"        : "${tumor_only.tumor_bam}"
    def targets_file     = target_bed     ? "--target ${target_bed}"       : ""
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""


    """
    freebayes \\
        -f $fasta \\
        $targets_file \\
        $samples_file \\
        $populations_file \\
        $cnv_file \\
        $args \\
        $input > ${prefix}.vcf

    bgzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
