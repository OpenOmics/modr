# Quantification rules
rule nanocount:
    """
    Data-processing step to get transcript counts.
    This tool is specifically built for ONT direct-RNA 
    data and it is supported by MultiQC.
    Github: https://github.com/a-slide/NanoCount
    @Input:
        Sorted Transcriptomic BAM file (scatter)
    @Output:
        NanoCount counts file
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

rule nanocount_aggregate:
    """
    Data-processing step to aggreagte per-sample raw and normalized 
    counts into a counts matrix.
    Github: https://github.com/a-slide/NanoCount
    @Input:
        NanoCount counts file (gather)
    @Output:
        NanoCount estimated raw counts matrix,
        NanoCount TPM normalized counts matrix 
    """
    input:
        counts = expand(join(workpath, "{name}", "counts", "{name}.nanocount.transcripts.tsv"), name=samples),
    output:
        est_count = join(workpath, "project", "counts", "known", "nanocount.transcripts.counts.tsv"),
        tpm = join(workpath, "project", "counts", "known", "nanocount.transcripts.tpm.tsv"),
    params:
        rname   = "nanocaggr",
        script  = join("workflow", "scripts", "create_matrix.py"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "nanocount_aggregate", cluster))
    shell: """
    # Estimated counts matrix,
    # raw counts 
    {params.script} \\
        --input {input.counts} \\
        --output {output.est_count} \\
        --join-on transcript_name \\
        --extract est_count \\
        --clean-suffix '.nanocount.transcripts.tsv' \\
        --nan-values 0.0
    # TPM counts matrix,
    # normalized counts
    {params.script} \\
        --input {input.counts} \\
        --output {output.tpm} \\
        --join-on transcript_name \\
        --extract tpm \\
        --clean-suffix '.nanocount.transcripts.tsv' \\
        --nan-values 0.0
    """


rule flair_correct:
    """
    Data-processing step to convert BAM file to BED12 format and 
    to correct misaligned splice sites using genome annotations.
    We are skipping the `flair align` step as we already have 
    a sorted BAM file aligned to the genome using the same 
    minimap2 options.
    Github: https://github.com/BrooksLabUCSC/flair
    @Input:
        Sorted Genomic BAM file (scatter)
    @Output:
        Sorted Genomics Alignments in BED12,
        FLAIR Correct Genomic Alignments in BED12  
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam"),
        ref = join(workpath, "refs", ref_genome),
        gtf = join(workpath, "refs", ref_gtf),
    output:
        bed12     = join(workpath, "{name}", "bams", "{name}.sorted.genome.bed"),
        corrected = join(workpath, "{name}", "bams", "{name}_all_corrected.bed"),
    params:
        rname  = "flaircorr",
        prefix = join(workpath, "{name}", "bams", "{name}"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['flair'], use_singularity)
    threads: int(allocated("threads", "flair_correct", cluster))
    shell: """
    # Convert BAM into BED12
    bam2Bed12 \\
        -i {input.bam} \\
    > {output.bed12}
    # Correct misaligned splice
    # sites with flair correct 
    flair correct \\
        -q {output.bed12} \\
        -g {input.ref} \\
        -f {input.gtf} \\
        -t {threads} \\
        -o {params.prefix} \\
        --nvrna
    """