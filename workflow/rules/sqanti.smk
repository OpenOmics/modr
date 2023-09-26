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
        --filter_mono_exonic \\
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


rule sqanti_nanocount:
    """
    Data-processing step to get known and novel sqanti-aligned transcript counts.
    For more info, please checkout: https://github.com/a-slide/NanoCount
    @Input:
        Sorted SQANTI-aligned Transcriptomic BAM file (scatter)
    @Output:
        NanoCount SQANTI-aligned (known/novel) counts file
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.sqanti.transcriptome.bam"),
    output:
        counts = join(workpath, "{name}", "counts", "{name}.nanocount.sqanti_filtered.transcripts.tsv"),
    params:
        rname   = "sqanti_nanocount",
        em_iter = config['options']['nanocount_em_iter'],
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "sqanti_nanocount", cluster))
    shell: """
    NanoCount \\
        -i {input.bam} \\
        --extra_tx_info \\
        -e {params.em_iter} \\
        -o {output.counts}
    """



rule sqanti_nanocount_aggregate:
    """
    Data-processing step to aggreagte per-sample raw and normalized 
    counts into a counts matrix. This counts matrix was generated via
    alignment to a SQANTI-filter derived transcriptome. As so, it will
    contain known/novel transcripts.
    Github: https://github.com/a-slide/NanoCount
    @Input:
        NanoCount SQANTI-aligned (known/novel) counts file (gather)
    @Output:
        NanoCount SQANTI-aligned (known/novel) estimated raw counts matrix,
        NanoCount SQANTI-aligned (known/novel) TPM normalized counts matrix 
    """
    input:
        counts = expand(join(workpath, "{name}", "counts", "{name}.nanocount.sqanti_filtered.transcripts.tsv"), name=samples),
    output:
        est_count = join(workpath, "project", "counts", "novel", "nanocount.sqanti_filtered.transcripts.counts.tsv"),
        tpm = join(workpath, "project", "counts", "novel", "nanocount.sqanti_filtered.transcripts.tpm.tsv"),
    params:
        rname   = "sqnati_nanocaggr",
        script  = join("workflow", "scripts", "create_matrix.py"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "sqanti_nanocount_aggregate", cluster))
    shell: """
    # Estimated counts matrix,
    # raw counts 
    {params.script} \\
        --input {input.counts} \\
        --output {output.est_count} \\
        --join-on transcript_name \\
        --extract est_count \\
        --clean-suffix '.nanocount.sqanti_filtered.transcripts.tsv' \\
        --nan-values 0.0
    # TPM counts matrix,
    # normalized counts
    {params.script} \\
        --input {input.counts} \\
        --output {output.tpm} \\
        --join-on transcript_name \\
        --extract tpm \\
        --clean-suffix '.nanocount.sqanti_filtered.transcripts.tsv' \\
        --nan-values 0.0
    """
