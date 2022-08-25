process RANDOMSITES {
    //tag "$meta.id"
    label 'process_low'

    input:
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
    path "versions.yml"                   , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "randomsites"

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
            -s ${RANDOM} \\ //SET RANDOM FUNCTION FROM BASH
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            snv > "random_snv_${prefix}.bed"
        """
    }

    else if (type == 'indel'){
        """
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -s ${RANDOM} \\ //SET RANDOM FUNCTION FROM BASH
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            indel --maxlen $maxlen > "random_indels_${prefix}.bed"
        """
    }

    else if (type == 'both'){
        """
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -s ${RANDOM} \\ //SET RANDOM FUNCTION FROM BASH
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            snv > "random_snv${prefix}.bed"

        
        python3 ${bamsurgeon_path}/scripts/randomsites.py \\
            $args \\
            -g $fasta \\
            -b $bed \\
            -s ${RANDOM} \\ //SET RANDOM FUNCTION FROM BASH
            -n $mut_number \\
            --minvaf $minvaf \\
            --maxvaf $maxvaf \\
            indel --maxlen $maxlen > "random_indels_${prefix}.bed"
        """
    }

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'random_sites.py'
        Type of variants generated: $type
        Number of variants generated: $mut_number
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: $bed
        FASTA used: $fasta
    END_VERSIONS
    """

    stub: //CHECK IF THIS SECTION IS MANDATORY
    def prefix = task.ext.prefix ?: "neat"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    END_VERSIONS
    """
}
