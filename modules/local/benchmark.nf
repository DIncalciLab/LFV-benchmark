process BENCHMARK {
    tag "Benchmark of the spiked-in somatic variants on sample: ${meta.sample}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf_vardict), path(vcf_mutect), path(vcf_varscan)

    output:
    tuple val(meta), path("*.png")   , emit: benchmark
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    python3 benchmark.py \\
        -n $neat \\
        -g $bamsurgeon \\
        -v $vcf_vardict \\
        -m $vcf_mutect \\
        -s $vcf_varscan

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        
    END_VERSIONS
    """
}