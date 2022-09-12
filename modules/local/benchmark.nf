process BENCHMARK {

    tag "Benchmark of the spiked-in somatic variants on sample: ${meta[0]}"
    label 'process_low'

    container "aldosr/neat:3.2"

    input:
    val(meta)

    output:
    path("*.png")   , emit: benchmark
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    benchmark.py \\
        -n "${meta[1]}" \\
        -g  \\
        -v  \\
        -m  \\
        -s 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        
    END_VERSIONS
    """
}