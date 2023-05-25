# Differential alternative splicing analysis rules
rule flair_altsplice:
    """
    Data-processing step to find differential alternative 
    splicing events. The `flair diffSplice` module requires
    at least 3 replicates in each group. As so, this rule 
    should not run if one group has 2 replicates. If a group
    has less than 3 replicates, the `diffsplice_fishers_exact` 
    script should be run instead. For more information, please
    read through flair's documenation:
    https://flair.readthedocs.io/en/latest/modules.html#flair-diffsplice
    Github: https://github.com/BrooksLabUCSC/flair
    @Input:
        Flair collapse isoforms BED file
        Flair transcript counts (scatter-per-contrast)
    @Output:
        Quantified alternative 3' splice sites (TSV),
        Quantified alternative 5' splice sites (TSV),
        Quantified skipped exon splice sites (TSV),
        Quantified retained intron splice sites (TSV),
        Differential alternative 3' splice sites (TSV),
        Differential alternative 5' splice sites (TSV),
        Differential skipped exon splice sites (TSV),
        Differential retained intron splice sites (TSV),
    """
    input:
        isoforms = join(workpath, "project", "counts", "novel", "flair.isoforms.bed"),
        counts   = join(workpath, "project", "counts", "novel", "flair.transcripts.counts.tsv"),
    output:
        a3   = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "diffsplice.alt3.events.quant.tsv"),
        a5   = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "diffsplice.alt5.events.quant.tsv"),
        es   = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "diffsplice.es.events.quant.tsv"),
        ir   = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "diffsplice.ir.events.quant.tsv"),
        dfa3 = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "drimseq_alt3_{grp1}_v_{grp2}.tsv"),
        dfa5 = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "drimseq_alt5_{grp1}_v_{grp2}.tsv"),
        dfes = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "drimseq_es_{grp1}_v_{grp2}.tsv"),
        dfir = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}", "drimseq_ir_{grp1}_v_{grp2}.tsv"),
    params:
        rname  = "flairaltsplice",
        outdir = join(workpath, "project", "diffsplice", "{grp1}_v_{grp2}"),
        group1 = "{grp1}",
        group2 = "{grp2}",
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['flair'], use_singularity)
    threads: int(allocated("threads", "flair_altsplice", cluster))
    shell: """
    # Find differential alternative
    # splicing events such as skipped
    # exons, alternative 5'/3' splice 
    # sites, and retained introns 
    flair diffSplice \\
        --threads {threads} \\
        --test \\
        --isoforms {input.isoforms} \\
        --counts_matrix {input.counts} \\
        --conditionA {params.group1} \\
        --conditionB {params.group2} \\
        --out_dir {params.outdir} \\
        --out_dir_force
    """