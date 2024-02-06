# Quality-control rules
# Pre-alignment quality control
rule pycoqc_summary:
    """
    Data-processing step to produce a summary file for pycoQC.
    The summary file is generated from a directory containing 
    basecalled fast5 files.
    pycoQC Github: https://github.com/a-slide/pycoQC
    @Output:
        Summary text file that is needed to run pycoQC. 
    """
    output:
        tsv  = join(workpath, "reports", "pycoQC", "{samples_hash}", "summary_sequencing.tsv"),
    params:
        rname  = "pyqcsum", 
        f5path = lambda w: hash_to_fast5[w.samples_hash]
    conda: depending(join(workpath, config['conda']['pycoqc']), use_conda)
    container: depending(config['images']['pycoqc'], use_singularity)
    threads: int(allocated("threads", "pycoqc_summary", cluster))
    shell: """
    # Create pycoQC summary file
    Fast5_to_seq_summary \\
        --threads {threads} \\
        -f {params.f5path} \\
        -s {output.tsv} \\
        --verbose 2 \\
        --fields barcode_arrangement barcode_full_arrangement \\
            barcode_score calibration_strand_end \\
            calibration_strand_genome_template \\
            calibration_strand_identity calibration_strand_start \\
            called_events channel channel_digitisation \\
            channel_offset channel_range channel_sampling_rate \\
            device_id duration flow_cell_id mean_qscore_template \\
            protocol_run_id read_id read_number run_id sample_id \\
            sequence_length_template skip_prob start_mux \\
            start_time stay_prob step_prob strand_score
    """


rule pycoqc_report:
    """
    Reporting step to produce a pycoQC report. pycoQC computes 
    several ONT-specific QC metrics and generates interactive QC 
    plots from a sequencing summary report.
    pycoQC Github: https://github.com/a-slide/pycoQC
    @Input:
        Summary text file that is needed to run pycoQC. 
    @Output:
        pycoQC HTML report
    """
    input:
        tsv  = join(workpath, "reports", "pycoQC", "{samples_hash}", "summary_sequencing.tsv"),
    output:
        html = join(workpath, "reports", "pycoQC", "{samples_hash}", "pycoQC_report.html"),
        json = join(workpath, "reports", "pycoQC", "{samples_hash}", "pycoQC_report.json"),
    params:
        rname  = "pyqcrpt", 
        f5path = lambda w: hash_to_fast5[w.samples_hash]
    conda: depending(join(workpath, config['conda']['pycoqc']), use_conda)
    container: depending(config['images']['pycoqc'], use_singularity)
    threads: int(allocated("threads", "pycoqc_report", cluster))
    shell: """
    pycoQC \\
        -f {input.tsv} \\
        -o {output.html} \\
        -j {output.json} 
    """


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
        tmpdir = tmpdir,
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "fastqc_raw", cluster))
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff
    fastqc \\
        -t {threads} \\
        -o "${{tmp}}" \\
        {input.fq}
    
    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\; 
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
        tmpdir = tmpdir,
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "fastqc_filtered", cluster))
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff
    fastqc \\
        -t {threads} \\
        -o "${{tmp}}" \\
        {input.fq}

    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\; 
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
        -o {params.outdir} \\
        --no_static \\
        --verbose
    """


rule nanostat:
    """
    Quality-control step to gather various statistics from a BAM file.
    This tool is supported by MultiQC.
    Github: https://github.com/wdecoster/nanostat
    @Input:
        Sorted Genomic BAM file (scatter)
    @Output:
        NanoStat metrics file
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam"),
    output:
        metrics = join(workpath, "{name}", "bams", "{name}.sorted.genome.metrics"),
    params:
        rname  = "nanostat",
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "nanostat", cluster))
    shell: """
    NanoStat \\
        -t {threads} \\
        --bam {input.bam} \\
    > {output.metrics}
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
        expand(join(workpath, "{name}", "bams", "{name}.sorted.genome.metrics"), name=samples),
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
