process BAMSURGEON_RANDOMSITES {
    tag "Create artificial random mutations for ${meta.sample}"
    label 'process_low'

    input:
    val meta
    val mut_number
    val minvaf
    val maxvaf
    val type
    val maxlen
    path fasta
    path bed
    path bamsurgeon_path

    output:
    tuple val(meta), path("*.bed")        , emit: bed
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

    if (type == 'snv') {
        """
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            snv > "${prefix}_snv"
        """
    } else if (type == 'indel'){
        """
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            indel --maxlen $maxlen > "${prefix}_indels"
        """
    } else {
        """
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            snv > "${prefix}_snv"

        
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            indel --maxlen $maxlen > "${prefix}_indels"
        """
    }

    """
    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py'
        Type of variants generated: $type
        Number of variants generated: $mut_number
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
    END_VERSIONS
    """
}