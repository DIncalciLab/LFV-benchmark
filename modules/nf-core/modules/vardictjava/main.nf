process VARDICTJAVA {
    tag "Variant calling using VarDict on BAMSurgeon spiked-in sample: ${meta.sample}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::vardict-java=1.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.2--hdfd78af_3':
        'quay.io/biocontainers/vardict-java:1.8.2--hdfd78af_3' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    
    path   bed
    path   fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf_vardict
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "vardict"
    def VERSION = '1.8.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.


    """
    vardict-java \
        $args \
        -b $bam \
        -th $task.cpus \
        -N $prefix \
        -G $fasta \
        $bed \
        | teststrandbias.R \
        | var2vcf_valid.pl \
            $args2 \
            -N $prefix \
        | gzip -c > ${prefix}_vardict.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: $VERSION
        var2vcf_valid.pl: \$(echo \$(var2vcf_valid.pl -h | sed -n 2p | awk '{ print \$2 }'))
    END_VERSIONS
    """

}