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
        manifest = join(workpath, "project", "counts", "novel", "flair.sample_manifest.tsv"),
    output:
        counts = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "flair.transcripts.counts.tsv"),
        deg = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "workdir", "genes_deseq2_{grp1}_v_{grp2}_results.tsv"),
        dei = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "workdir", "isoforms_deseq2_{grp1}_v_{grp2}_results.tsv"),
        diu = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}", "workdir", "isoforms_drimseq_{grp1}_v_{grp2}_results.tsv"),
    params:
        rname   = "flairdiffexp",
        outdir  = join(workpath, "project", "diffexp", "{grp1}_v_{grp2}"),
        group1  = "{grp1}",
        group2  = "{grp2}",
        script  = join("workflow", "scripts", "cut2tsv.py"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['flair'], use_singularity)
    threads: int(allocated("threads", "flair_diffexp", cluster))
    shell: """
    # Subset counts matrix to contain 
    # only samples belonging to the 
    # contrast of interest, this will
    # ensure the correct group is set
    # as the baseline for design matrix.
    # By default, the first group it 
    # encounters in the sample names 
    # will be the intercept. 
    # Get the first group's samples
    grp1_samples=$(
        awk -F '\\t' -v OFS='_' \\
            '$2=="{params.group1}" {{print $1,$2,$3}}' \\
            {input.manifest}
    )
    # Get the second group's samples
    grp2_samples=$(
        awk -F '\\t' -v OFS='_' \\
            '$2=="{params.group2}" {{print $1,$2,$3}}' \\
            {input.manifest}
    )
    # Create per-contrast counts matrix,
    # extracts samples by their col name
    {params.script} -c ids ${{grp1_samples}} ${{grp2_samples}} \\
        -i {input.counts} \\
    > {output.counts}

    # Find differential gene/isoform 
    # expression and find differential
    # isoform usage
    flair diffExp \\
        --threads {threads} \\
        --counts_matrix {output.counts} \\
        --out_dir {params.outdir} \\
        --exp_thresh 10 \\
        --out_dir_force
    """
