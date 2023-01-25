# Quantification rules
rule nanocount:
    """
    Quality-control step to gather various statistics from a BAM file.
    This tool is supported by MultiQC.
    Github: https://github.com/a-slide/NanoCount
    @Input:
        Sorted Transcriptomic BAM file (scatter)
    @Output:
        NanoCount metrics file
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.transcriptome.bam"),
    output:
        counts = join(workpath, "{name}", "counts", "{name}.nanocount.transcripts.tsv"),
    params:
        rname   = "nanocount",
        em_iter = config['options']['nanocount_em_iter'],
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "nanocount", cluster))
    shell: """
    NanoCount \\
        -i {input.bam} \\
        --extra_tx_info \\
        -e {params.em_iter} \\
        -o {output.counts} 
    """
