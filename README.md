## Introduction

**dincalcilab/lowfrac-variant-benchmark** is a bioinformatics pipeline to generate syntethic data sets to benchmark low-fraction somatic variant callers.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have installed from the [nf-core/modules](https://github.com/nf-core/modules) repository.

## Pipeline summary

1. Generate artificial BAM files with [`NEAT`](https://github.com/ncsa/NEAT)
2. Generate random variant sites (SNP and indels) with [`BAMsurgeon`](https://github.com/adamewing/bamsurgeon)
3. Spike-in the random variants in the artificial BAM files and generate a ground truth VCF ([`BAMsurgeon`](https://github.com/adamewing/bamsurgeon))
4. Benchmark three variant callers against the ground truth:
  - [`MuTect 2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
  - [`VarDict`](https://github.com/AstraZeneca-NGS/VarDictJava)
  - [`VarScan2`](https://dkoboldt.github.io/varscan/)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)
2. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run dincalcilab/lowfrac-variant-benchmark -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   Example:

   ```console
   nextflow run dincalcilab/lowfrac-variant-benchmark --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Credits

dincalcilab/lowfrac-variant-benchmark was originally written by Aldo Sergi and Luca Beltrame.
## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  dincalcilab/lowfrac-variant-benchmark for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
