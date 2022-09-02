process BAMSURGEON {
    
    tag "Spike-in artificial random mutain in sample: ${meta.sample}"
    label 'process_high'
   
    container "aldosr/bamsurgeon:1.3"
    
    input:
    val meta
    val mut_number
    val minvaf
    val maxvaf
    val maxlen
    val fasta
    val type

    path bed
    path fasta
    path picardjar

    output:
    tuple val(meta), path("*.txt")        , emit: random_mut
    tuple val(meta), path("*.yml")        , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
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
    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g "${fasta}" \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        snv > ${prefix}_random_snv.txt

    
    python3 -O bamsurgeon/bin/addsnv.py \
        $args2 \
        -v "random_mut/${prefix}_random_snv.txt" \
        -f "${meta.info}" \
        -r "${fasta}" \
        -o "spiked_snv/${prefix}_spiked_snv.bam" \
        --picardjar $picardjar \
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp_addsnv" 


    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_snv.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Only SNVs were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: "${bed}"
        FASTA used: "${fasta}"
        Meta: $meta
    END_VERSIONS
    """

    } else if (type == 'indel'){

    """
    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g "${fasta}" \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        indel --maxlen $maxlen > "random_mut/${prefix}_random_indel.txt"


    python 3 -O bamsurgeon/bin/addindel.py \
        $args3 \
        -v "random_mut/${prefix}_random_indel.txt" \
        -f "${meta.info}" \
        -r "${fasta}" \
        -o "spiked_indel/${prefix}_spiked_indel.bam" \
        --picardjar $picardjar \
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp_addindel"

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_indel.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Only INDELs were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        Max INDELs length: $maxlen
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS
    """

    } else if (type == 'both') {

    """
    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g "${fasta}" \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        snv > "random_mut/${prefix}_random_snv.txt"

    python3 bamsurgeon/scripts/randomsites.py \
        $args \
        -g "${fasta}" \
        -b $bed \
        -n $mut_number \
        --minvaf $minvaf \
        --maxvaf $maxvaf \
        indel --maxlen $maxlen > "random_mut/${prefix}_random_indel.txt"

    
    python3 -O bamsurgeon/bin/addsnv.py \
        $args2 \
        -v "random_mut/${prefix}_random_snv.txt" \
        -f ${meta.info} \ 
        -r "${fasta}" \
        -o "spiked_snv/${prefix}_spiked_snv.bam" \
        --picardjar $picardjar \
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp_addsnv"
    
    python 3 -O bamsurgeon/bin/addindel.py \
        $args3 \
        -v "random_mut/${prefix}_random_indel.txt" \
        -f "${meta.info}" \
        -r "${fasta}" \
        -o "spiked_snv_indel/${prefix}_spiked_snv_indel.bam" \
        --picardjar $picardjar \
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp_addindel"


    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_snv.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Both SNVs and INDELs (Max Length = $maxlen) were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: "${bed}"
        FASTA used: "${fasta}"
        Meta: $meta
    END_VERSIONS
    """
    
    } else {
        log.info 'ERROR: YOU MUST SPECIFY A MUTATION TYPE TO SPIKEIN'
    }

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
    END_VERSIONS
    """
}
