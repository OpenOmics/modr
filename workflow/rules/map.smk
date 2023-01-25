# Data processing rules
rule minimap2_genome:
    """
    Data-processing step to align direct RNA reads against reference genome,
    also runs samtools stats and flagstat to get alignment statistics.
    @Input:
        Nanofilt quality filtered FastQ file (scatter),
        Genomic FASTA file
    @Output:
        Genomic BAM file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
        ref = join(workpath, "refs", ref_genome),
    output:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam"),
        bai = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam.bai"),
        flagstat = join(workpath, "{name}", "bams", "{name}.sorted.genome.flagstats"),
        stats    = join(workpath, "{name}", "bams", "{name}.sorted.genome.stats"),
        idx      = join(workpath, "{name}", "bams", "{name}.sorted.genome.idxstats"),
    params:
        rname  = 'minimap2',
        tmpdir = join(workpath, "{name}", "bams", "genome_tmp"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "minimap2_genome", cluster)) 
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
        # conversion of U to T bps 
        minimap2 \\
            -ax splice \\
            -uf \\
            -k14 \\
            {input.ref} \\
            {input.fq} \\
        | samtools sort -@{threads} \\
            -T "${{tmp}}" \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            -
        
        # Gather BAM statistics
        samtools idxstats {output.bam} > {output.idx}
        samtools flagstat {output.bam} > {output.flagstat}
        samtools stats {output.bam} > {output.stats}
        """


rule minimap2_transcriptome:
    """
    Data-processing step to align direct RNA reads against transcriptome.
    @Input:
        Nanofilt quality filtered FastQ file (scatter),
        Transcriptomic FASTA file
    @Output:
        Transcriptomic BAM file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.filtered.fastq.gz"),
        ref = join(workpath, "refs", ref_transcripts),
    output:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.transcriptome.bam"),
        bai = join(workpath, "{name}", "bams", "{name}.sorted.transcriptome.bam.bai"),
    params:
        rname = 'minimap2',
        tmpdir = join(workpath, "{name}", "bams", "transcriptome_tmp"),
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "minimap2_transcriptome", cluster)) 
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
        # conversion of U to T bps
        minimap2 \\
            -ax splice \\
            -N 10 \\
            -uf \\
            -k14 \\
            {input.ref} \\
            {input.fq} \\
        | samtools sort -@{threads} \\
            -T "${{tmp}}" \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            -
        """


rule stranded_bigwigs:
    """
    Data-processing step to create stranded bigWigs tracks for IGV.
    Using deeptools: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
    @Input:
        Genomic BAM file (scatter),
    @Output:
        Normalized (TPM) BigWig of fwd strand,
        Normalized (TPM) BigWig of rev strand 
    """
    input:
        bam = join(workpath, "{name}", "bams", "{name}.sorted.genome.bam"),
    output:
        fwd = join(workpath, "{name}", "bigwigs", "{name}.tpm.fwd.bw"),
        rev = join(workpath, "{name}", "bigwigs", "{name}.tpm.rev.bw"),
    params:
        rname = 'bigwigs',
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "stranded_bigwigs", cluster)) 
    shell: 
        """
        # TPM normalized BigWig
        # of the forward strand
        bamCoverage \\
            -b {input.bam} \\
            -o {output.fwd} \\
            --samFlagExclude 16 \\
            --normalizeUsing BPM

        # TPM normalized BigWig
        # of the reverse strand
        bamCoverage \\
            -b {input.bam} \\
            -o {output.rev} \\
            --samFlagInclude 16 \\
            --normalizeUsing BPM
        """