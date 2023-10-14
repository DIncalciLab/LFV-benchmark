## Introduction

**dincalcilab/LFV-benchmark** is a bioinformatics pipeline to generate syntethic data sets to benchmark low-fraction somatic variant callers (in case, variants with higher fractions can be benchmarked).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple computing infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes should be installed from the [nf-core/modules](https://github.com/nf-core/modules) repository.

## Pipeline summary

1. Generate artificial BAM files with [`NEAT`](https://github.com/ncsa/NEAT)
2. Generate random variant sites (SNV and indels) with [`BAMsurgeon`](https://github.com/adamewing/bamsurgeon)
3. Spike-in the random variants in the artificial BAM files and generate a ground truth VCF ([`BAMsurgeon`](https://github.com/adamewing/bamsurgeon))
4. Benchmark the following variant callers against the ground truth:
  - [`MuTect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  - [`VarDict`](https://github.com/AstraZeneca-NGS/VarDictJava)
  - [`VarScan2`](https://dkoboldt.github.io/varscan/)
  - [`LoFreq`](https://github.com/CSB5/lofreq)
  - [`FreeBayes`](https://github.com/freebayes/freebayes)
  - [`Strelka2`](https://github.com/Illumina/strelka)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)
2. Install [`Picard`](https://github.com/broadinstitute/picard), mandatory for running BAMsurgeon
3. Download the pipeline:

   ```console
   nextflow pull https://github.com/DIncalciLab/LFV-benchmark
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Test the pipeline. Different use cases are given below.


## Use cases

### Generate 2 tumor/normal pair with a coverage of 100x with random SNV/INDEL
  
   ```console
   nextflow run DIncalciLab/LFV-benchmark --outdir <OUTDIR> --fasta <FASTA> --picardjar <PICARDJAR> --samples 2 --coverage 100
   ```
  setting the path for the output folder to OUTDIR, the path of the fasta file (e.g. hg19/38) to FASTA and the path of the picard jar file (generally in `build/libs/picard.jar`, see the [`Picard repo`](https://github.com/broadinstitute/picard)) to PICARDJAR.

### Perform the benchmark on N tumor/normal pair given in input to the pipeline:
    
   ```console
   nextflow run DIncalciLab/LFV-benchmark --outdir <OUTDIR> --fasta <FASTA> --picardjar <PICARDJAR> --input_normal <NORMALBAM> --input_tumor <TUMORBAM> --skip_normal_generation --skip_tumor_generation
   ```
   and give to `--input_normal` and `--input_tumor` the path of the folder with the test BAM files, normal and tumor respectively (can be found in the `test` folder of this repo).
   
## Usage

The pipeline steps can be run either together or separately. To run the entire workflow (e.g. generate artificial datasets, spike-in somatic variants and benchmark the variant callers), run:

```console
nextflow run DIncalciLab/LFV-benchmark \
        --outdir OUTDIR \
        [--readlen INT] \
        --fasta FASTA \
        [--bed BED] \
        [--coverage INT] \
        [--error_model ERROR_MODEL] \
        [--mutation_model MUT_MODEL] \
        [--gc_model GC_MODEL] \
        [--fraglen_model FRAGLEN_MODEL] \
        --picardjar PICARDJAR \
        [--type snv] \
        --samples 5

mandatory arguments:
  --outdir              output directory
  --fasta               path to the fasta file
  --picardjar           path to the picard.jar file
  --samples             nunber of artificial samples to be generated
  
optional arguments:
  --readlen             int (default: 151)        length of artificial reads generated by NEAT
  --bed                 FASTA                     print version and exit
  --coverage            int (default: 10000)      Coverage depth of the artificial samples generated by NEAT
  --error_model         ERROR_MODEL               path of sequencing error model used by NEAT
  --mutation_model      MUT_MODEL                 path of mutational rate model used by NEAT
  --gc_model            GC_MODEL                  path of GC model used by NEAT
  --fraglen_model       FRAGLEN_MODEL             path of the fragment length model used by NEAT
  --type                string (default: both)    type of low-fraction variants to spike-in in the artificial samples. Can be: snv, indel, both         
```
Default models utilized by NEAT are available at: https://github.com/ncsa/NEAT/tree/master/models

### Skip the generation of normal samples

To spike-in variants in existing normal/tumor samples, just skip the NEAT step (i.e. the generation of artificial normal samples) adding the following commands to the pipeline:

```console
--input_normal              path to the folder containing the BAM samples to spike-in and benchmark
--skip_normal_generation
```

### Skip the generation of normal and tumor samples

To run only the variant calling benchmark on existing normal/tumor BAM files, the spike-in step (operated by BAMsurgeon) can be skipped, using the following command:

```console
[--input_normal]          path to the folder containing the normal BAM samples to be used in the tumor-normal paired mode 
--input_tumor             path to the folder containing the BAM samples to benchmark
--skip_normal_generation
--skip_tumor_generation
```

### Skip the generation of tumor samples

To run the benchmark on the artificial files generated by NEAT (for testing purpouse mainly), skip only the spike-in process using the following command:

```console
--skip_tumor_generation
```

### Skip the benchmark step

The benchmark step can be skipped using the following command:

```console
--skip_benchmark
```

Moreover, each variant caller can be skipped independently, using the commands `--skip_vardict`, `--skip_mutect` and so on...

### High-sensitivity mode

To run the variant callers using the optimized set of parameters (see the related paper for details), just add the command `--high-sensitivity` to the pipeline.

## Credits

DIncalciLab/LFV-benchmark was originally written by Aldo Sergi and Luca Beltrame.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  DIncalciLab/LFV-benchmark for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
