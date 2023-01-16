process STRELKA_SOMATIC {
    tag 'Variant calling using Strelka2 (somatic) on BAMSurgeon spiked-in sample: ${meta.sample_name}'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta), path(normal),   path(tumor)
    path(manta_candidate_small_indels)//, path(manta_candidate_small_indels_tbi),
          //path(target_bed),              path(target_bed_index)
    val fasta
    val target_bed
    //path  fai

    output:
    //tuple val(meta),
    path("*.somatic_indels.vcf.gz")    , emit: vcf_strelka_indels
    //tuple val(meta),
    path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_strelka_indels_tbi
    //tuple val(meta),
    path("*.somatic_snvs.vcf.gz")      , emit: vcf_strelka_snvs
    //tuple val(meta),
    path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_strelka_snvs_tbi
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_name}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : "" //disabled due to bug in strelka2
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""

    if ( params.high_sensitivity ){
    """
    awk  "BEGIN{FS=OFS="="}/maxIndelSize/{${2}=90}1" /usr/local/share/strelka-2.9.10-1/bin/configureStrelkaSomaticWorkflow.py.ini

    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumor.tumor_bam} \\
        --normalBam ${normal.normal_bam} \\
        --referenceFasta $fasta \\
        ${options_manta} \\
        $args \\
        --runDir strelka
    """
    }
    if ( !params.high_sensitivity ){
    """
    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumor.tumor_bam} \\
        --normalBam ${normal.normal_bam} \\
        --referenceFasta $fasta \\
        ${options_manta} \\
        $args \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j $task.cpus

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
    }
}
