process GENERATE_PLOTS {

    tag "Generate plots and calculate benchmark of the spiked-in somatic variants"
    label 'process_low'

    container "aldosr/cyvcf2:0.30.18"

    input:
    path(neat)
    path(bamsurgeon)
    val(vardict)
    val(mutect)
    val(varscan)

    output:
    path("*.xlsx")                   , optional: true, emit: benchmark
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ""
    
    """
    benchmark.py \\
        -n ${neat} \\
        -b ${bamsurgeon} \\
        -v ${vardict} \\
        -m ${mutect} \\
        -s ${varscan}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        
    END_VERSIONS
    """
}