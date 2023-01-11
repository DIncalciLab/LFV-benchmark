process VARDICTJAVA {
    tag "Variant calling using VarDict on BAMSurgeon spiked-in sample: ${meta.sample_name}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::vardict-java=1.8.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.2--hdfd78af_3':
        'quay.io/biocontainers/vardict-java:1.8.2--hdfd78af_3' }"

    input:
    tuple val(meta), path(normal),  path(tumor)
    
    val   fasta
    path  bed

    output:
    //tuple val(meta),
    path("*.vcf.gz"),                  emit: vcf_vardict
    //path("*.vcf.gz.tbi"),              emit: vcf_vardict_tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "vardict"
    def bam     = ( normal_bam  ) ? "'${tumor_bam}|${normal_bam}'" : "${tumor_only.tumor_bam}"
    def mode = ( !(${normal_bam} == 'EMPTY')  ) ?
                "testsomatic.R | var2vcf_paired.pl -N '${prefix}_tumor|${prefix}_normal'" :
                "teststrandbias.R | var2vcf_valid.pl -N ${prefix} -E"
    def VERSION = '1.8.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def avail_mem = 3
    if (!task.memory) {
        log.info '[VarDict Java] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }



   if ( params.high_sensitivity ){
   """
   vardict-java \
       -G ${fasta} \
       -f 0.0001 \
       -N ${prefix} \
       -b $bam \
       $args \
       $bed \
           | $mode  -f 0.0001 | gzip -c > ${prefix}.vcf.gz

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       vardict-java: $VERSION
       var2vcf_valid.pl: \$(echo \$(var2vcf_valid.pl -h | sed -n 2p | awk '{ print \$2 }'))
       var2vcf_paired.pl: \$(echo \$(var2vcf_paired.pl -h | sed -n 2p | awk '{ print \$2 }'))
   END_VERSIONS
   """
   }

   if ( !params.high_sensitivity ) {
   """
   vardict-java \
       -G ${fasta} \
       -f 0.01 \
       -N ${prefix} \
       -b $bam \
       $args \
       $bed \
           | $mode  -f 0.01 | gzip -c > ${prefix}.vcf.gz

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       vardict-java: $VERSION
       var2vcf_valid.pl: \$(echo \$(var2vcf_valid.pl -h | sed -n 2p | awk '{ print \$2 }'))
       var2vcf_paired.pl: \$(echo \$(var2vcf_paired.pl -h | sed -n 2p | awk '{ print \$2 }'))
   END_VERSIONS
   """
   }

}