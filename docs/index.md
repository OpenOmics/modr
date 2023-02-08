<div align="center">

  <h1 style="font-size: 250%">modr 🔬</h1>

  <b><i>RNA modification Pipeline for ONT direct RNA-sequencing</i></b><br>
  <a href="https://hub.docker.com/repository/docker/skchronicles/dinopore">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/dinopore">
  </a>
  <a href="https://github.com/OpenOmics/modr/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/modr/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/modr/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/modr/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/modr/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/modr?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/modr/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/modr">
  </a>

  <p>
    This is the home of the pipeline, modr. Its long-term goals: to detect RNA modifications in Oxford Nanopore direct RNA-sequencing data like no pipeline before!
  </p>

</div>  

## Overview
Welcome to modr's documentation! This guide is the main source of documentation for users that are getting started with the [ONT RNA modification pipeline](https://github.com/OpenOmics/modr/). 

The **`./modr`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">modr <b>run</b></code>](usage/run.md)   
    Run the modr pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">modr <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">modr <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">modr <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**modr** is a pipeline to detect RNA modifications in Oxford Nanopore direct RNA sequencing data. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with direct RNA-sequencing data generated from [Oxford Nanopore sequencing Technologies](https://nanoporetech.com/). As input, it accepts a set of gzipped FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/modr/issues).

## Contribute 

This site is a living document, created for and by members like you. modr is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/modr).

<!---
## Citation

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    Coming soon!
    ```

=== "APA"

    ```
    Coming soon!
    ```

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.7477631).
--->

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
