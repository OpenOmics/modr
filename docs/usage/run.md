# <code>modr <b>run</b></code>

## 1. About 
The `modr` executable is composed of several inter-related sub commands. Please see `modr -h` for all available options.

This part of the documentation describes options and concepts for <code>modr <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running modr pipeline. 

Setting up the modr pipeline is fast and easy! In its most basic form, <code>modr <b>run</b></code> only has *four required inputs*.

## 2. Synopsis
```text
$ modr run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {slurm,local}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--resource-bundle RESOURCE_BUNDLE] [--use-conda] \
      [--contrasts CONTRASTS] [--quality-filter QUALITY_FILTER] \
      [--rna-editing {A-I, ...}] \
      --genome {hg38_41, mm10_M25, mm39_M31} \
      --groups GROUPS \
      --input INPUT [INPUT ...] \
      --output OUTPUT
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, a reference genome via the `--genome` argument, and a groups file containing additional sample metadata via the `--groups` argument.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input Oxford Nanopore FastQ files(s).**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed. If a sample has multiple fastq files for different barcodes, the pipeline expects each barcoded FastQ file endwith the following extension: `_N.fastq.gz`, where `N` is a number. Internally, the pipeline will concatenate each of these FastQ files prior to processing the data. Here is an example of an input sample with multiple barcode sequences: `S1_0.fastq.gz`, `S1_1.fastq.gz`, `S1_2.fastq.gz`, `S1_3.fastq.gz`. Given this barcoded sample, the pipeline will create the following concatenated FastQ file: `S1.fastq.gz`. 
> 
> ***Example:*** `--input .tests/*.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/modr_out`

---  
  `--genome {hg38_41, mm10_M25, mm39_M31}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. modr does comes bundled with pre-built reference files from GENCODE for human and mouse samples. Please select from one of the following options: `hg38_41`, `mm10_M25`, `mm39_M31`. Please note that `hg38_41` is a human reference genome, while `mm10_M25` and `mm39_M31` are two different reference genomes available for mouse.
> ***Example:*** `Example: --genome hg38_41`

---  
  `--groups GROUPS`
> **Groups file.**   
> *type: file*
>   
> This tab-delimited (TSV) file is used to collect additional metadata that is required to run the pipeline. This file consists of three columns containing  the names of each sample, the path to the sample's fast5 directory, and the name of a group a sample is assoicated with. Each of these fields is required! The first column should contain the base name of a given sample. The base name of a given sample can be determined by removing its file extension from the sample's FastQ file, for example: `WT_S4.fastq.gz` becomes `WT_S4` in the groups file. Here is an example base name of a multiplexed sample. 
> 
> Given:
> `WT_S2_1.fastq.gz`, `WT_S2_2.fastq.gz`, `WT_S2_3.fastq.gz`
> The base name would be `WT_S2`.
>
> The second column should contain the path to the sample's fast5 directory. The fast5 directory contains HDF5 files which couple sequencing information with raw event signal from the ONT sequencer. The raw signal information is used to infer RNA modifcation. As so, the pipeline needs access to both the FastQ files and their fast5 files. The third, and last column, contains group information for each sample. A user can set their own groups. A group can represent anything from a timepoint, to a experimental condition, to a treatment, etc. Please note that groups are composed of alphanumeric characters, cannot startwith a number, and cannot contain any `-` characters. Groups can contain `_` character.
>
> *Contents of example groups file:*
> ```
> WT_S1	.tests/WT_S1/fast5/	T1
> WT_S2	.tests/WT_S2/fast5/	T1
> WT_S3	.tests/WT_S3/fast5/	T2
> WT_S4	.tests/WT_S4/fast5/	T2
> WT_S5	.tests/WT_S5/fast5/	T2
> ```

> ***Example:*** `Example: --groups .tests/groups.tsv`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 


  `--contrasts CONTRASTS`
> **Contrasts file for Differential Analyses.**   
> *type: file*
>   
> This tab-delimited (TSV) file is used to setup comparisons between different sets of groups. These comparsions are used to find differential transcript expression (DTE), differential transcript usage (DTU), and differential alternative splicing (AS) events. Please see the `--groups` option above for more information about how to define groups within a set of samples. The contrasts file consists of two columns containing the names of each group to compare. The names defined in this file must also exist in the groups file.
> 
> Given the following groups file:
> ```
> WT_S1	.tests/WT_S1/fast5/	T1
> WT_S2	.tests/WT_S2/fast5/	T1
> WT_S3	.tests/WT_S3/fast5/	T2
> WT_S4	.tests/WT_S4/fast5/	T2
> WT_S5	.tests/WT_S5/fast5/	T2
> ```
>
> Here is an example contrasts file:
> ```
> T2	T1
> ```
>
> **Please note:** the order of the groups defined in the contrasts file will determine how to interpret the direction of the fold-change. As so, if we use the `T2	 T1` comparison as an example for interpreting DTE results, a positive fold-change for a given transcript would indicate that the `T2` samples' expression is higher than the `T1` samples' expression. 

> ***Example:*** `Example: --contrasts .tests/contrasts.tsv`

--- 
  `--quality-filter QUALITY_FILTER`  
> **Quality score filter.**  
> *type: int*
> *default: 8*
> 
> This option filters reads on a minimum average quality score. Any reads with an average minimum quality score less than this threshold will be removed. The default average minimum quality filter is set to 8.
> 
> ***Example:*** `--quality-filter 10`

--- 
  `--rna-editing {A-I, ...}`  
> **Type of RNA editing or RNA modifications to identify.**  
> *type: str*
> *default: None*
> 
> RNA editing has been show to affect the function and stability of RNA molecules. This option allow a user to define the type of RNA editing sites to identify. 
>
> The pipeline can currently identify the following types of RNA editing (more coming soon!):  

> - **A-I**: A common type of RNA editing catalysed by double-stranded RNA (dsRNA)-specific adenosine deaminase (ADAR) enzymes.  

> 
> ***Example:*** `--rna-editing A-I`

### 2.3 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running modr in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:modr*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:modr".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `modr cache` subcommand can be used to create a local SIF cache. Please see `modr cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running modr with this option when ever possible.
> 
> ***Example:*** `--singularity-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`

---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOBID`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /scratch/$USER/`

---  
  `--resource-bundle RESOURCE_BUNDLE`
> **Path to a resource bundle downloaded with the install sub command.**  
> *type: path*  
>
> The resource bundle contains the set of required reference files for processing any data. The path provided to this option will be the path to the `modr` directory that was created when running the install sub command. Please see the install sub command for more information about downloading the pipeline's resource bundle.
> 
> ***Example:*** `--resource-bundle /data/$USER/refs/modr`

---  
  `--use-conda`   
> **Use Conda/mamba instead of Singularity.**  
> *type: boolean flag*
> 
> Use Conda/Mamba instead of Singularity. By default, the pipeline uses singularity for handling required software dependencies. This option overrides that behavior, and it will use Conda or mamba instead of Singularity. The use of Singuarity and Conda are mutually exclusive. Please note that conda or mamba must be in your $PATH prior to running the pipeline.
> 
> ***Example:*** `--use-conda`

### 2.4 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./modr run --input .tests/*.fastq.gz \
           --output test_01 \
           --genome hg38_41 \
           --groups .tests/groups.tsv \
           --mode slurm \
           --dry-run

# Step 2B.) Run the modr pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./modr run --input .tests/*.fastq.gz \
           --output test_01 \
           --genome hg38_41 \
           --groups .tests/groups.tsv \
           --mode slurm
```
