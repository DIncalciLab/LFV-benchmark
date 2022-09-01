process BAMSURGEON_SPIKEIN {
    tag "Spike-in artificial mutation in ${meta.sample} sample"
    label 'process_high'

    container 'aldosr/bamsurgeon:1.3'

    input:
    //tuple val(meta) file(bam) UNCOMMENT WHEN USING NEAT
    file
    tuple val(meta) file(snv)
    tuple val(meta) file(snv)
    val type
    path fasta
    path picardjar

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.vcf")        , emit: vcf
    path "${meta.sample}_versions.yml"                   , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def version = '1.3' //VERSION IS HARDCODED

    def avail_mem = 3
    if (!task.memory) {
        log.info '[BAMSurgeon] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (type == 'snv') {
    """
    python3 -O bamsurgeon/bin/addsnv.py \
        $args \
        -v $snv \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp" 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/add_snv.py'
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS
    """
    }
    
    else if (type == 'indel') {
    """
    python 3 -O bamsurgeon/bin/addindel.py \\
        $args2 \\
        -v $indel \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/add_indel.py'
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS 
    """
    }

    else if (type == 'both') {
    """
    python3 -O bamsurgeon/bin/addsnv.py \\
        $args \\
        -v $snv \
        -f $bam \
        -r $fasta \
        -o "add_snv/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp"

        python3 -O bamsurgeon/bin/addindel.py \\
        $args \\
        -v $indel \
        -f $bam \
        -r $fasta \
        -o "add_indel/snv_${prefix}.bam" \
        --picardjar $picardjar \\ 
        --alignopts c:250,M:,t:$task.cpus,v:1 \
        -p $task.cpus \
        --tmpdir "tmp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/add_snv.py' + 'BAMSurgeon/add_indel.py'
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS
    """
    }
}
