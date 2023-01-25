# Quality-control rules
# Pre-alignment quality control
rule fastqc_raw:
    """
    Quality-control step to assess sequencing quality of raw data.
    @Input:
        Setup FastQ file (scatter)
    @Output:
        FastQC html report and zip file on raw data
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.fastq.gz"),
    output:
        html     = join(workpath, "{name}", "fastqc", "{name}_fastqc.html"),
        zip_file = join(workpath, "{name}", "fastqc", "{name}_fastqc.zip")
    params:
        rname  = "rawfqc", 
        outdir = join(workpath, "{name}", "fastqc"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "fastqc_raw", cluster))
    shell: """
    fastqc \\
        -t {threads} \\
        -o {params.outdir} \\
        {input.fq}
    """


rule fastqc_filtered:
    """
    Quality-control step to assess sequencing quality of filtered data.
    @Input:
        Nanofilt filtered FastQ file (scatter)
    @Output:
        FastQC html report and zip file on filtered data
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
    output:
        html     = join(workpath, "{name}", "fastqc", "{name}.filtered_fastqc.html"),
        zip_file = join(workpath, "{name}", "fastqc", "{name}.filtered_fastqc.zip")
    params:
        rname  = "filtfqc", 
        outdir = join(workpath, "{name}", "fastqc"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "fastqc_filtered", cluster))
    shell: """
    fastqc \\
        -t {threads} \\
        -o {params.outdir} \\
        {input.fq}
    """


# Post-alignment quality control
rule nanoplot:
    """
    Quality-control step to visualize QC information.
    Github: https://github.com/wdecoster/NanoPlot
    @Input:
        Sorted Genomic BAM file (scatter)
    @Output:
        Nanoplot html report and various QC plots
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam"),
    output:
        html = join(workpath, "{name}", "nanoplot", "NanoPlot-report.html"),
    params:
        rname  = "nanoplot", 
        outdir = join(workpath, "{name}", "nanoplot"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "nanoplot", cluster))
    shell: """
    NanoPlot \\
        -t {threads} \\
        --bam {input.bam} \\
        -o {params.outdir}
    """


rule multiqc:
    """
    Reporting step to aggregate QC information across
    all supported tools.
    Docs: https://multiqc.info/
    @Input:
        FastQC html reports (gather)

    @Output:
        MulitQC html report
    """
    input:
        expand(join(workpath, "{name}", "fastqc", "{name}_fastqc.zip"), name=samples),
        expand(join(workpath, "{name}", "fastqc", "{name}.filtered_fastqc.zip"), name=samples),
        expand(join(workpath, "{name}", "bams", "{name}.sorted.genome.stats"), name=samples),
    output:
        html = join(workpath, "reports", "multiqc_report.html"),
    params:
        rname  = "multiqc", 
        wdir   = join(workpath),
        outdir = join(workpath, "reports"),
        conf   = join(workpath, "resources", "multiqc_config.yaml"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "mulitqc", cluster))
    shell: """
    multiqc \\
        --ignore '*/.singularity/*' \\
        -f \\
        -c {params.conf} \\
        --interactive \\
        --outdir {params.outdir} \\
        {params.wdir}
    """