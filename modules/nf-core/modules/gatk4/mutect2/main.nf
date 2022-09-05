process GATK4_MUTECT2 {
    tag "Variant calling using VarDict on BAMSurgeon spiked-in sample: ${meta.sample}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4-4.2.3.0-1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_1':
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    
    val fasta
    path bed

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf_mutect
    tuple val(meta), path("*.tbi")        , emit: tbi_mutect
    tuple val(meta), path("*.stats")      , emit: stats_mutect
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "mutect2"
    //def inputs = input.collect{ "--input $it"}.join(" ")
    //def interval_command = intervals ? "--intervals $intervals" : ""
    //def pon_command = panel_of_normals ? "--panel-of-normals $panel_of_normals" : ""
    //def gr_command = germline_resource ? "--germline-resource $germline_resource" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (mode == 'high-sensitivity'){
    """
    gatk --java-options "-Xmx${avail_mem}g" Mutect2 \\
        -I $bam \\
        --reference $fasta \\
        -L $bed \\
        -O ${prefix}.vcf \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    } else{
    """
        gatk --java-options "-Xmx${avail_mem}g" Mutect2 \\
        -I $bam \\
        --reference $fasta \\
        -L $bed \\
        $args \\
        -O ${prefix}.vcf \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

    """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
