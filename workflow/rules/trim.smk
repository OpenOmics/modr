# Helper functions
def get_barcoded_fastqs(wildcards):
    """
    Returns a list of per-sample multiplexed barcoded FastQ files.
    """
    barcodes = samples2barcodes[wildcards.name]
    # Convert strings to integers,
    # prior to sorting on keys,
    # JSON cannot have int keys.
    # {'0': 'WT_S1_0.fastq.gz', '1': 'WT_S1_1.fastq.gz'} 
    # { 0:  'WT_S1_0.fastq.gz',  1: 'WT_S1_1.fastq.gz'}
    barcodes = {int(k):v for k,v in barcodes.items()}
    if barcodes:
        # Merge multiple barcoded FastQ files,
        # Sort files based on barcode int to
        # ensure output is deterministic
        sorted_keys = sorted(barcodes.keys())              # 0, 1, ...
        sorted_values = [barcodes[k] for k in sorted_keys] # WT_S1_0.fastq.gz, WT_S1_1.fastq.gz
        return [join(workpath, f) for f in sorted_values]
    else:
        # Already merged, return input file
        return ["{0}.fastq.gz".format(join(workpath, wildcards.name))]


# Data processing rules
rule setup:
    """
    Initialization step to either demultiplex samples with multiple
    barcodes, or create symlinks to already merged inputs. Th shell 
    command is dynamically resolved depened on whether a given sample
    has a set of barcode files (meaning it needs concatenation) or if
    the sample is already merged (no barcode files, create symlink). 
    @Input:
        FastQ file (gather-per-sample-multiple-barcodes)
    @Output:
        Merged/symlinked FastQ files
    """
    input:
        fq=get_barcoded_fastqs,
    output:
        fq=join(workpath, "{name}", "fastqs", "{name}.fastq.gz"),
    params:
        rname ='merge',
         # Building to merge multiple files or 
         # create a symlink if already merged
        prefix = lambda w: "cat" 
            if samples2barcodes[w.name] 
            else "ln -rsf",
        suffix = lambda w: ">" 
            if samples2barcodes[w.name] 
            else "",
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    shell: 
        """
        {params.prefix} {input.fq} {params.suffix} {output.fq}
        """


rule nanofilt:
    """
    Data-processing step to perform base quality filtering with Nanofilt. 
    @Input:
        Setup FastQ file (scatter)
    @Output:
        Quality filtered FastQ file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.fastq.gz"),
    output:
        flt = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
    params:
        rname='nanofilt',
        qual_filt=quality_filter,
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    shell: 
        """
        # Nanofilt requires uncompressed input
        gunzip -c {input.fq} \\
            | NanoFilt -q {params.qual_filt} \\
            | gzip \\
        > {output.flt}
        """
