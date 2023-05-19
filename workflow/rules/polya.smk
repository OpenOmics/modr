# Poly-A tail length estimation rules
rule nanopolish_polya:
    """
    Data-processing step to get poly-A tail lengths
    from nanopolish. This tool is specifically built 
    for ONT direct-RNA data. For more information,
    please checkout it's documentation:
    https://nanopolish.readthedocs.io/en/latest/quickstart_polya.html
    Github: https://github.com/jts/nanopolish
    @Input:
        Filtered FasstQ file (scatter)
        Sorted Transcriptomic BAM file (scatter)
        Transcriptome fasta file 
    @Output:
        Nanopolish polyA tail length estimates
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
        bam = join(workpath, "{name}", "bams", "{name}.sorted.transcriptome.bam"),
        ref = join(workpath, "refs", ref_transcripts),
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