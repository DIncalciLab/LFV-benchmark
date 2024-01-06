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
2. Install [`Picard`](https://github.com/broadinstitute/picard) (`<3.00` to avoid errors with Java 17), mandatory for running BAMsurgeon
3. Download the pipeline:

   ```console
   nextflow pull https://github.com/DIncalciLab/LFV-benchmark
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`. A simple configuration (`test_local`) to run the pipeline in local on a low-memory machine is also provided. You can change it depending on your needs.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Test the pipeline. Different use cases are given below.
5. NB: Due to a bug in [`Pandas/Numpy old versions `]([https://nf-co.re/tools/#downloading-pipelines-for-offline-use](https://github.com/pandas-dev/pandas/issues/39520)) it is advised to run the script to calculate variant calling performance separately from the pipeline, as described in the "Use cases" section.


## Use cases

### Generate 2 tumor/normal pair with a coverage of 30x with random SNV

   ```console
   nextflow run DIncalciLab/LFV-benchmark --input_all <input_csv> --outdir <OUTDIR> --bed <path_to_bed> --fasta <FASTA> --picardjar <PICARDJAR> --samples 2 --coverage 30 --type snv --high-sensitivity -profile test_local
   ```
  where:
```console
  --input_all           path to csv file containing the name of the samples to be generated (`/assets/samplesheet.csv` is an example of this file)
  --outdir              path of the output folder
  --bed                 path of the bed file, containing the regions to be generated
  --fasta               path of the fasta file (e.g. hg19/hg38). NB: the fasta folder should aldo contain the fasta index file
  --picard              path of the picard.jar file (generally in `build/libs/picard.jar`, see the [`Picard repo`](https://github.com/broadinstitute/picard))
  --samples             number of samples to be generated
  --coverage            coverage of the artificial samples to be generated
  --type                type of mutations spiked-in (SNVs/INDELs/both)
  --high-sensitivity    use tuned parameters for variant calling
  ```
  Then calculate the performance of the callers by launching the script `benchmark_standalone.py` located in the `/bin` folder, using the following command:

   ```console
   ./benchmark_standalone.py -t snv -s <BAMSURGEON_VCF> -v <VARIANT_CALLING_FOLDER> --o <OUTDIR>
   ```
  where:
```console
  --t              type of variants inserted from BAMSURGEON
  --s              directory with VCF files generated from BAMSURGEON with the random SNVs inserted in the generated normal samples
  --v              output folder generated from the pipeline with the variant calling outputs (normally `outdir/variant_calling`
  --o              output folder
```


### Perform the benchmark on 2 existent tumor/normal pair with a coverage of 100X and with random SNVs inserted:

   ```console
   nextflow run DIncalciLab/LFV-benchmark --outdir <OUTDIR> --fasta <FASTA> --picardjar <PICARDJAR> --input_normal <NORMALBAM> --input_tumor <TUMORBAM> --skip_normal_generation --skip_tumor_generation --high-sensitivity
   ```
   and give to `--input_normal` and `--input_tumor` the path of the folder with the test BAM files, normal and tumor respectively (can be found in the `test_files` folder of this repo).

   Then calculate the performance as described previously. The VCF files generated from BAMSURGEON with the inserted SNVs can be found in `test_files/spiked_vcf/SNV/100X`

### Perform the benchmark on a single tumor sample with a coverage of 30.000X with spiked SNVs/INDELs from BAMSURGEON
   ```console
   nextflow run DIncalciLab/LFV-benchmark --outdir <OUTDIR> --fasta <FASTA> --picardjar <PICARDJAR> --input_tumor <TUMORBAM> --skip_normal_generation --skip_tumor_generation --high-sensitivity --tumor_only
   ```
Input tumor can be found in `test_files/tumor_bam/high_coverage`
Then calculate the performance of the callers by launching the script `benchmark_standalone_tumor_only.py` located in the `/bin` folder, using the following command:

   ```console
   ./benchmark_standalone_tumor_only.py -t both -s <BAMSURGEON_VCF> -i <BAMSURGEON_VCF_INDEL> -v <VARIANT_CALLING_FOLDER> --o <OUTDIR>
   ```
   where `--i` is the directory with VCF files generated from BAMSURGEON with the random INDELs inserted in the generated samples. The VCF files generated from BAMSURGEON with the inserted SNVs/INDELs can be found in `test_files/spiked_vcf/SNV/high_coverage` and `test_files/spiked_vcf/INDEL/high_coverage`

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
