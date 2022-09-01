process BAMSURGEON_RANDOMSITES {
    tag "Create artificial random mutations for ${meta.sample}"
    label 'process_low'

   //conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
   /*
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'aldosr/bamsurgeon:1.3' }"*/
    
    input:
    val meta
    val mut_number
    val minvaf
    val maxvaf
    val maxlen
    path fasta
    path bed
    path bamsurgeon_path

    output:
    tuple val(meta), path("${prefix}.txt")        , emit: mut
    tuple val(meta), path("*.yml")        , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def version = '1.3' //VERSION IS HARDCODED
    
    def avail_mem = 3
    if (!task.memory) {
        log.info '[BAMSurgeon/random_sites.py] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    touch ${prefix}.txt
    """
/*
    """
    echo "test" > test.txt
    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g $fasta \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        snv > "${prefix}_snv.txt"

    
    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g $fasta \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        indel --maxlen $maxlen > "${prefix}_indel.txt"
    """
*/
    """
    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py'
        Number of variants generated: $mut_number
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: $bed
        FASTA used: $fasta
        Meta: $prefix
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
    END_VERSIONS
    """
}
