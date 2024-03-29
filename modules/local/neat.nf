process NEAT {
    tag "Create artificial normal datasets for sample: ${meta.sample_name}"
    label 'process_high'

    container "aldosr/neat:3.2"

    input:
    val meta
    path bed
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz")                                   , emit: vcf

    tuple val(meta), path("*.bam"), path("*.bai")                       , emit: bam

    tuple val(meta), path("*.yml")                                      , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def error_model = params.error_model ?: ''
    def mut_model = params.mutation_model ?: ''
    def gc_model = params.gc_model ?: ''
    def fraglen_model = params.fraglen_model ?: ''
    def prefix = task.ext.prefix ?: ''
    def version = '3.2' //VERSION IS HARDCODED


    def avail_mem = 3
    if (!task.memory) {
        log.info '[NEAT] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gen_reads.py \
        $args \
        -r $fasta \
        -R ${params.readlen} \
        --pe-model $fraglen_model \
        -c ${params.coverage} \
        -e $error_model \
        --gc-model $gc_model \
        -tr $bed \
        -m $mut_model \
        -o $prefix

    samtools index ${prefix}_golden.bam

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        ncsa/NEAT: 'Version $version'
        GC bias model: $gc_model
        Bed used: $bed
        Mutational model: $mut_model
        FASTA: $fasta
        Sequencing error model: $error_model
        Frag length model: $fraglen_model
        Coverage: ${params.coverage}
        Read length: ${params.readlen}
        Meta: $meta
    END_VERSIONS
     """

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.bam

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        ncsa/NEAT: 'Version $version'
        GC bias model: $gc_model
        Bed used: $bed
        Mutational model: $mut_model
        FASTA: $fasta
        Sequencing error model: $error_model
        Frag length model: $fraglen_model
        Coverage: ${params.coverage}
        Read length: ${params.readlen}
    END_VERSIONS
    """
}
