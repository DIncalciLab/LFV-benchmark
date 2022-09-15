process GENERATE_PLOTS {

    tag "Generate plots and calculate benchmark of the spiked-in somatic variants"
    label 'process_low'

    container "aldosr/neat:3.2"

    input:
    val(vardict)

    output:
    path("*.xlsx")                   , emit: benchmark
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    benchmark.py \\
        -n ${meta[0].vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        
    END_VERSIONS
    """
}