process BAMSURGEON {
    //tag "$meta.id"
    label 'process_high'

    input:
    path bam
    path bed
    path fasta
    path bamsurgeon_path
    path picardjar

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.vcf")        , emit: vcf
    path "versions.yml"                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "bamsurgeon"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[BAMSurgeon] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (spike_type == 'snv') {
    """
    python3 -O ${bamsurgeon_path}/bin/addsnv.py \\
        $args \\
        -v $bed \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp" 
    """
    }
    
    else if (spike_type == 'indel') {
    """
    python 3 -O ${bamsurgeon_path}/bin/addindel.py \\
        $args \\
        -v $bed \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp" 
    """
    }

    else if (spike_type == 'both') {
    """
    python3 -O ${bamsurgeon_path}/bin/addsnv.py \\
        $args \\
        -v $bed \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp"

        python3 -O ${bamsurgeon_path}/bin/addindel.py \\
        $args \\
        -v $bed \
        -f $bam \
        -r $fasta \
        -o "add_indel/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp"
    """
    }

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon'
        BED used: $bed
        FASTA used: $fasta
    END_VERSIONS
    """    
}
