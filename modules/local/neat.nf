process NEAT {
    tag "$meta.id"
    label 'process_high'

    input:
    val readlen
    val coverage
    path fasta
    path fraglenmodel
    path seqerrormodel
    path gcbiasmodel
    path bed
    path mutmodel

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.bam")        , emit: bam
    path "versions.yml"                   , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = '3.2' //VERSION IS HARDCODED


    def avail_mem = 3
    if (!task.memory) {
        log.info '[ncsa/NEAT] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }


    """
    python3 /path/to/NEAT/gen_reads.py \\
        $args \\
        -r $fasta \\
        -R $readlen \\
        --pe-model $fraglenmodel \\
        -c $coverage \\
        -e $seqerrormodel \\
        --gc-model $gcbiasmodel \\
        -tr $bed \\
        --rng ${RANDOM} \\ //COME PASSARE RANDOM DI BASH?
        -m $mutmodel \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncsa/NEAT: 'Version $version'
        GC bias model: $gcbiasmodel
        Bed used: $bed
        Mutational model: $mutmodel
        Mutational rate: none   
        FASTA: $fasta
        Sequencing error model: $seqerrormodel
        Frag length model: $fraglenmodel
        Coverage: $coverage
        Read length: $readlen
    END_VERSIONS
     """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncsa/NEAT: 'Version $version'
        GC bias model: $gcbiasmodel
        Bed used: $bed
        Mutational model: $mutmodel
        Mutational rate: none   
        FASTA: $fasta
        Sequencing error model: $seqerrormodel
        Frag length model: $fraglenmodel
        Coverage: $coverage
        Read length: $readlen
    END_VERSIONS
    """
    }
