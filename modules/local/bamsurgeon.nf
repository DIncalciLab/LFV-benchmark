process BAMSURGEON {

    tag "Spike-in artificial somatic mutations in sample: ${meta.sample_name}"
    label 'process_medium'

    container "aldosr/bamsurgeon:1.3-custom"

    input:
    tuple val(meta), val(neat)

    val (mut_number)
    val (minvaf)
    val (maxvaf)
    val (maxlen)
    val (fasta)

    path (bed)
    path (picardjar)

    output:
    tuple val(meta), path("*.txt")                                                 , emit: random_mut

    tuple val(meta), path("bamsurgeon*.bam"), path("bamsurgeon*.bai")              , emit: bam

    tuple val(meta), path("*.vcf.gz"), path("*.tbi")                               , emit: vcf

    tuple val(meta), path("*.yml")                                                 , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: ''
    def version = '1.3' //VERSION IS HARDCODED

    def avail_mem = 3
    if (!task.memory) {
        log.info '[BAMSurgeon/random_sites.py] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    if (params.type == 'snv') {

    """
    randomsites.py \\
        $args \\
        -g "${fasta}" \\
        -b $bed \\
        -n $mut_number \\
        --minvaf $minvaf \\
        --maxvaf $maxvaf \\
        snv > ${prefix}_random_snv.txt


    addsnv.py \\
        $args2 \\
        -v ${prefix}_random_snv.txt \\
        -f "${neat}" \\
        -r "${fasta}" \\
        -o ${prefix}_spiked_snv.bam \\
        --picardjar $picardjar \\
        --alignopts c:250,M:,t:$task.cpus,v:1 \\
        -p $task.cpus \\
        --tmpdir tmp_addsnv

    samtools index ${prefix}_spiked_snv.bam

    bcftools reheader --fai "${fasta}.fai" ${prefix}_spiked_snv.addsnv.${prefix}_random_snv.vcf \
        | bcftools sort | bgzip -c > ${prefix}_spiked_snv.vcf.gz
    tabix -p vcf ${prefix}_spiked_snv.vcf.gz

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_snv.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Only SNVs were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: "${bed}"
        FASTA used: "${fasta}"
        Meta: $meta
    END_VERSIONS
    """

    } else if (params.type == 'indel'){

    """
    randomsites.py \\
        $args \\
        -g "${fasta}" \\
        -b $bed \\
        -n $mut_number \\
        --minvaf $minvaf \\
        --maxvaf $maxvaf \\
        indel --maxlen $maxlen > ${prefix}_random_indel.txt


    addindel.py \\
        $args3 \\
        -v ${prefix}_random_indel.txt \\
        -f "${neat}" \\
        -r "${fasta}" \\
        -o ${prefix}_spiked_indel.bam \\
        --picardjar $picardjar \\
        --alignopts c:250,M:,t:$task.cpus,v:1 \\
        -p $task.cpus \\
        --tmpdir tmp_addindel

    samtools index ${prefix}_spiked_indel.bam

    bcftools reheader --fai "${fasta}.fai" ${prefix}_spiked_indel.addindel.${prefix}_random_indel.vcf \
        | bcftools sort | bgzip -c > ${prefix}_spiked_indel.vcf.gz
    tabix -p vcf ${prefix}_spiked_indel.vcf.gz

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_indel.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Only INDELs were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        Max INDELs length: $maxlen
        BED used: $bed
        FASTA used: $fasta
        Meta: $meta
    END_VERSIONS
    """

    } else if (params.type == 'both') {

    """
    randomsites.py \\
        $args \\
        -g "${fasta}" \\
        -b $bed \\
        -n $mut_number \\
        --minvaf $minvaf \\
        --maxvaf $maxvaf \\
        snv > ${prefix}_random_snv.txt

    randomsites.py \\
        $args \\
        -g "${fasta}" \\
        -b $bed \\
        -n $mut_number \\
        --minvaf $minvaf \\
        --maxvaf $maxvaf \\
        indel --maxlen $maxlen > ${prefix}_random_indel.txt


   addsnv.py \\
        $args2 \\
        -v ${prefix}_random_snv.txt \\
        -f "${neat}" \\
        -r "${fasta}" \\
        -o ${prefix}_spiked_snv.bam \\
        --picardjar $picardjar \\
        --alignopts c:250,M:,t:$task.cpus,v:1 \\
        -p $task.cpus \\
        --tmpdir tmp_addsnv

    samtools index ${prefix}_spiked_snv.bam

    addindel.py \\
        $args3 \\
        -v ${prefix}_random_indel.txt \\
        -f ${prefix}_spiked_snv.bam \\
        -r "${fasta}" \\
        -o ${prefix}_spiked_snv_indel.prereplace \\
        --picardjar $picardjar \\
        --alignopts c:250,M:,t:$task.cpus,v:1 \\
        -p $task.cpus \\
        --tmpdir tmp_addindel

    java -jar $picardjar AddOrReplaceReadGroups I=${prefix}_spiked_snv_indel.prereplace \
        O=${prefix}_spiked_snv_indel.bam \\
        VALIDATION_STRINGENCY=LENIENT \\
        RGID=NEAT \\
        RGLB=NEAT \\
        RGPL=NEAT \\
        RGPU=NEAT \\
        RGSM=NEAT

    samtools index ${prefix}_spiked_snv_indel.bam

    bcftools reheader --fai "${fasta}.fai" ${prefix}_spiked_snv_indel.addindel.${prefix}_random_indel.vcf \
        | bcftools sort | bgzip -c > ${prefix}_spiked_snv_indel.vcf.gz
    tabix -p vcf ${prefix}_spiked_snv_indel.vcf.gz

    mv ${prefix}_spiked_snv.bam ${prefix}_spiked_snv.bam.snv
    mv ${prefix}_spiked_snv.bam.bai ${prefix}_spiked_snv.bam.bai.snv

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
        BAMSurgeon: 'Version $version'
        Script: 'BAMSurgeon/random_sites.py + BAMSurgeon/add_snv.py'
        Number of variants generated: $mut_number
        Type of mutations spiked-in: 'Both SNVs and INDELs (Max Length = $maxlen) were inserted'
        Min VAF: $minvaf
        Max VAF: $maxvaf
        BED used: "${bed}"
        FASTA used: "${fasta}"
        Meta: $meta
    END_VERSIONS
    """

    } else {
        log.info 'ERROR: YOU MUST SPECIFY A MUTATION TYPE TO SPIKEIN'
    }

    stub:
    def prefix = task.ext.prefix ?: ""
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > "${prefix}.versions.yml"
    "${task.process}":
    END_VERSIONS
    """
}
