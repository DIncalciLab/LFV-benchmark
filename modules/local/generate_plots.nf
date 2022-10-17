process GENERATE_PLOTS {

    tag "Generate plots and calculate benchmark of the spiked-in somatic variants"
    label 'process_low'

    container "aldosr/cyvcf2:0.30.18"

    input:
    path(neat)
    path(bamsurgeon)
    path(vardict)
    path(mutect)
    path(varscan)

    output:
    path("*_variants.txt")                   , emit: variants_txt
    //path("*_neat.xlsx")                  , emit: neat_xlsx

    //path("*_bamsurgeon.txt")             , emit: bamsurgeon_txt
    //path("*_bamsurgeon.xlsx")            , emit: bamsurgeon_xlsx

    //path("*_vardict.txt")                , emit: vardict_txt
    //path("*_vardict.xlsx")               , emit: vardict_xlsx

    //path("*_mutect.txt")                 , emit: mutect_txt
    //path("*_mutect.xlsx")                , emit: mutect_xlsx

    //path("*_varscan.txt")                , emit: varscan_txt
    //ath("*_varscan.xlsx")               , emit: varscan_xlsx

    path("benchmark.txt")              ,  emit: benchmark_txt
    //path("*_benchmark.xlsx")             ,  emit: benchmark_xlsx

    path("versions.yml")                  , emit: versions

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