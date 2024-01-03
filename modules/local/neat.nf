process NEAT {
    tag "Create artificial normal datasets for sample: ${meta.sample_name}"
    label 'process_high'

    container "aldosr/neat:3.2"

    input:
    val meta
    val coverage
    path bed
    path fasta
    path neat_path

    output:
    tuple val(meta), path("*.vcf.gz")                    , emit: vcf
    //tuple val(meta), path("*.tbi")                       , emit: tbi

    tuple val(meta), path("*.bam")                       , emit: bam
    tuple val(meta), path("*.bai")                       , emit: bai

    tuple val(meta), path("*.yml")                       , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def error_model = params.error_model ? "-e $params.error_model" : ''
    def mut_model = params.mut_model ? "-m $params.mut_model" : ''
    def gc_model = params.gb_model ? "--gc-model $params.gc_model" : ''
    def fraglen_model = params.fraglen_model ? "--pe-model $params.fraglen_model" : ''
    def prefix = task.ext.prefix ?: ''
    def version = '3.2' //VERSION IS HARDCODED


    def avail_mem = 3
    if (!task.memory) {
        log.info '[NEAT] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    python3 ${neat_path}/gen_reads.py \
        $args \
        -r $fasta \
        -R $readlen \
        --pe-model $fraglenmodel \
        -c $coverage \
        -e $seqerrormodel \
        --gc-model $gcbiasmodel \
        -tr $bed \
        -m $mutmodel \
        -o $prefix

    samtools index ${prefix}_golden.bam

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        ncsa/NEAT: 'Version $version'
        GC bias model: $gcbiasmodel
        Bed used: $bed
        Mutational model: $mutmodel
        FASTA: $fasta
        Sequencing error model: $seqerrormodel
        Frag length model: $fraglenmodel
        Coverage: $coverage
        Read length: $readlen
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
        GC bias model: $gcbiasmodel
        Bed used: $bed
        Mutational model: $mutmodel
        FASTA: $fasta
        Sequencing error model: $seqerrormodel
        Frag length model: $fraglenmodel
        Coverage: $coverage
        Read length: $readlen
    END_VERSIONS
    """
}
