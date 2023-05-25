# Differential expression and differential usage rules
rule flair_diffexp:
    """
    Data-processing step to find differential expression 
    of genes/isoforms and to find differential isoform 
    usage. The `flair diffExp` module requires at least
    3 replicates in each group. This rule should not run 
    if one group has 2 replicates. If a group does have 
    less than 3 replicates, the `diff_iso_usage` script
    should be run instead. For more information, please
    read through flair's documenation:
    https://flair.readthedocs.io/en/latest/modules.html#flair-diffexp
    Github: https://github.com/BrooksLabUCSC/flair
    @Input:
        Flair transcript counts (scatter-per-contrast)
    @Output:
        Filtered differential gene expression table (TSV),
        Filtered differential isoform expression table (TSV),
        Filtered differential isoform usage table (TSV),
        Misc DESeq2 QC plots (PDF)
    """
    input:
        counts = join(workpath, "project", "counts", "novel", "flair.transcripts.counts.tsv"),
    output:
        deg = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "genes_deseq2_{grp1}_v_{grp2}.tsv"),
        dei = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "isoforms_deseq2_{grp1}_v_{grp2}.tsv"),
        diu = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "isoforms_drimseq_{grp1}_v_{grp2}.tsv"),
    params:
        rname   = "flairdiffexp",
        outdir  = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['flair'], use_singularity)
    threads: int(allocated("threads", "flair_diffexp", cluster))
    shell: """
    # Find differential gene/isoform 
    # expression and find differential
    # isoform usage
    flair diffExp \\
        --threads {threads} \\
        --counts_matrix {input.counts} \\
        --out_dir {params.outdir} \\
        --exp_thresh 10 \\
        --out_dir_force
    """