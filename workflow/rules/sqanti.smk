# Sqanti related quality-control and filtering rules,
# Sqanti is being used to annotate/characterize novel 
# isoforms and to build an even higher-confidence,
# filtered set of unique transcripts from flair.
# The resulting annotation/transcriptome will be
# used to for re-quantification of known AND novel
# isoforms using nanocount. Conceptually, this two-
# pass re-mapping approach is similar to STAR's. 
# In the first pass, we map reads against the ref
# genome, then we collapse/filter  to define high-
# confidence isoforms, then we take the resulting
# transcriptome for this project (containing high-
# confidence isoforms) and re-map it in the second
# pass.
rule sqanti_qc:
    """
    Data-processing step to characterize the input transcriptome 
    by computing a series of attributes by transcript, which are
    written to the classification file, and a series of attributes
    by junction, which are written to the junctions file. Please 
    note although we are running SQANTI3, the actual version of 
    the tool we are using is 'v5.1.2'. For more information, 
    please read through sqanti3's documenation:
    https://github.com/ConesaLab/SQANTI3/wiki/
    Github: https://github.com/ConesaLab/SQANTI3
    @Input:
        High-confidence Isoforms (FASTA) from flair collapse 
    @Output:
        Sqanti Classification file (TSV),
        Corrected Annotation (GTF),
        Corrected Transcriptome (FASTA)
    """
    input:
        fa     = join(workpath, "project", "counts", "novel", "flair.isoforms.fa"),
        gtf    = join(workpath, "refs", ref_gtf),
        genome = join(workpath, "refs", ref_genome),
    output:
        txt = join(workpath, "project", "counts", "novel", "sqanti.isoforms_classification.txt"),
        fa  = join(workpath, "project", "counts", "novel", "sqanti.isoforms_corrected.fasta"),
        gtf = join(workpath, "project", "counts", "novel", "sqanti.isoforms_corrected.gtf"),
    params:
        rname  = "sqanti_qc",
        prefix = "sqanti.isoforms",
        outdir      = join(workpath, "project", "counts", "novel"),
        cage_peak   = config['references'][genome]['SQANTI_CAGEPEAK'],
        polya_motif = config['references']['SQANTI_POLYA_MOTIF_LIST'],
    container: depending(config['images']['sqanti3'], use_singularity),
    threads: int(allocated("threads", "sqanti_qc", cluster)),
    shell: """
    # Structural and quality annotation
    # of collapsed transcript from flair
    sqanti3_qc.py \\
        {input.fa} \\
        {input.gtf} \\
        {input.genome} \\
        --aligner_choice=minimap2 \\
        --cpus {threads} \\
        --output {params.prefix} \\
        --dir {params.outdir} \\
        --CAGE_peak {params.cage_peak} \\
        --polyA_motif_list {params.polya_motif} \\
        --isoAnnotLite \\
        --force_id_ignore \\
        --fasta
    """


rule sqanti_ml_filter:
    """
    Data-processing step to filter the sqanti qc output. The author
    from sqanti highly recommends filtering its output before using
    it for down-stream analysis. Sqanti has a new filtering method
    that employs random forest to discriminate potential artifacts
    from true isoforms without the need for user-defined rules or
    manually-set thresholds (i.e. previous method). For more info, 
    please read through sqanti3's documenation:
    https://github.com/ConesaLab/SQANTI3/wiki/
    Github: https://github.com/ConesaLab/SQANTI3
    @Input:
        Sqanti Classification file (TSV),
        Corrected Annotation (GTF),
        Corrected Transcriptome (FASTA)
    @Input:
        ML Filtered Sqanti Classification file (TSV),
        ML Filtered Corrected Annotation (GTF),
        ML Filtered Corrected Transcriptome (FASTA)
    """
    input:
        txt = join(workpath, "project", "counts", "novel", "sqanti.isoforms_classification.txt"),
        fa  = join(workpath, "project", "counts", "novel", "sqanti.isoforms_corrected.fasta"),
        gtf = join(workpath, "project", "counts", "novel", "sqanti.isoforms_corrected.gtf"),
    output:
        txt = join(workpath, "project", "counts", "novel", "sqanti.isoforms_MLresult_classification.txt"),
        fa  = join(workpath, "project", "counts", "novel", "sqanti.isoforms.filtered.fasta"),
        gtf = join(workpath, "project", "counts", "novel", "sqanti.isoforms.filtered.gtf"),
    params:
        rname  = "sqanti_ml_filter",
        prefix = "sqanti.isoforms",
        outdir      = join(workpath, "project", "counts", "novel"),
    container: depending(config['images']['sqanti3'], use_singularity),
    threads: int(allocated("threads", "sqanti_ml_filter", cluster)),
    shell: """
    # Applies ML filter on selected
    # SQANTI3 QC classification file
    sqanti3_filter.py ML \\
        --output {params.prefix} \\
        --dir {params.outdir} \\
        --gtf {input.gtf} \\
        --isoforms {input.fa} \\
        {input.txt}
    """


rule sqanti_minimap2:
    """
    Data-processing step to align direct RNA reads against filtered, high-
    confidence annotated transcriptome from SQANTI.
    @Input:
        Nanofilt quality filtered FastQ file (scatter),
        Transcriptomic FASTA file
    @Output:
        Transcriptomic BAM file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
        ref  = join(workpath, "project", "counts", "novel", "sqanti.isoforms.filtered.fasta"),
    output:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.sqanti.transcriptome.bam"),
        bai = join(workpath, "{name}", "bams", "{name}.sorted.sqanti.transcriptome.bam.bai"),
    params:
        rname = 'sqanti_minimap2',
        tmpdir = join(workpath, "{name}", "bams", "sqanti_transcriptome_tmp"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "sqanti_minimap2", cluster)) 
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Align against reference genome,
        # minimap2 automatically handles
        # conversion of U to T bps, See 
        # this issue for the author of 
        # minimap2's recomendations for 
        # aligning direct RNA reads against 
        # the transcriptome:
        # https://github.com/lh3/minimap2/issues/340  
        minimap2 \\
            -ax map-ont \\
            -N 10 \\
            -k10 \\
            {input.ref} \\
            {input.fq} \\
        | samtools sort -@{threads} \\
            -T "${{tmp}}" \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            -
        """