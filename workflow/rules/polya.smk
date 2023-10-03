# Poly-A tail length estimation rules
rule nanopolish_index:
    """
    Data-processing step to build an index mapping 
    from base-called reads to the signals measured 
    by the sequencer. The `nanopolish polya` step 
    requires index files for the `--reads` option.
    For more information, please checkout its docs:
    https://nanopolish.readthedocs.io/en/latest/manual.html#index
    Github: https://github.com/jts/nanopolish
    @Input:
        Filtered FastQ file (scatter)
    @Output:
        Nanopolish polyA tail length estimates
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
    output:
        idx = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz.index"),
        rdb = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz.index.readdb"),
        gzi = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz.index.gzi"),
        fai = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz.index.fai"),
    params:
        rname   = "nanoindex",
        f5 = lambda w: str(config['fast5'][w.name]),
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "nanopolish_index", cluster))
    shell: """
    # Build an index to map basecalled
    # reads to ONT device's raw signal 
    nanopolish index \\
        --directory={params.f5} \\
        {input.fq}
    """


rule nanopolish_polya:
    """
    Data-processing step to get poly-A tail lengths
    from nanopolish. This tool is specifically built 
    for ONT direct-RNA data. For more information,
    please checkout its documentation:
    https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html
    Github: https://github.com/jts/nanopolish
    @Input:
        Filtered FastQ file (scatter)
        Sorted Transcriptomic BAM file (scatter)
        Transcriptome fasta file 
    @Output:
        Nanopolish polyA tail length estimates
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
        bam = join(workpath, "{name}", "bams", "{name}.sorted.transcriptome.bam"),
        ref = join(workpath, "refs", ref_transcripts),
        idx = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz.index"),
    output:
        polya = join(workpath, "{name}", "polyA", "{name}.nanopolish.transcripts.polyA.tsv"),
    params:
        rname   = "nanopolya",
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "nanopolish_polya", cluster))
    shell: """
    # Estimate polyA tail lengths
    # for each transcript
    nanopolish polya \\
        --threads={threads} \\
        --reads={input.fq} \\
        --bam={input.bam} \\
        --genome={input.ref} \\
    > {output.polya}
    """


rule nanopolish_polya_filter:
    """
    Data-processing step to filter and average the polyA
    tail length for each transcript. The estimated polyA 
    tail lengths need to be filtered on the qc_tag, only 
    keep estimates with the tag 'PASS'.
    @Input:
        Nanopolish polyA tail length estimates (scatter)
    @Output:
        Filtered/averaged nanopolish polyA tail length estimates
    """
    input:
        polya = join(workpath, "{name}", "polyA", "{name}.nanopolish.transcripts.polyA.tsv"),
    output:
        polya = join(workpath, "{name}", "polyA", "{name}.nanopolish.transcripts.polyA.average.tsv"),
    params:
        rname   = "filterpolya",
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "nanopolish_polya_filter", cluster))
    shell: """
    # Filter and average polyA tail length
    # of each transcript
    echo -e 'id\\tavg_polya_length' > {output.polya} 
    awk -F '\\t' -v OFS='\\t' 'NR>1 && $NF=="PASS" {{print}}' {input.polya} \\
        | cut -f2,9 \\
        | awk -F '\\t' -v OFS='\\t' \\
            '{{seen[$1]+=$2; count[$1]++}} END {{for (x in seen) print x, seen[x]/count[x]}}' \\
    >> {output.polya}
    """


rule nanopolish_polya_matrix:
    """
    Data-processing step to create a matrix containing
    the average polyA tail length for each transcript for
    each sample.
    @Input:
        Filtered/averaged nanopolish polyA tail length estimates (gather)
    @Output:
        Matrix containing average polyA tail lengths, 
        for each sample for each transcript (TSV)
    """
    input:
        polyas = expand(join(workpath, "{name}", "polyA", "{name}.nanopolish.transcripts.polyA.average.tsv"), name=samples),
    output:
        matrix = join(workpath, "project", "polyA", "nanopolish.transcripts.average_polyA_length.tsv"),
    params:
        rname   = "matrixpolya",
        script  = join("workflow", "scripts", "create_matrix.py"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "nanopolish_polya_matrix", cluster))
    shell: """
    # Create counts matrix each transcripts
    # average polyA tail length
    {params.script} \\
        --input {input.polyas} \\
        --output {output.matrix} \\
        --join-on id \\
        --extract avg_polya_length \\
        --clean-suffix '.nanopolish.transcripts.polyA.average.tsv' \\
        --nan-values 0.00
    """