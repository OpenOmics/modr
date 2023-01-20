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
    params:
        rname = 'minimap2',
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "minimap2_genome", cluster)) 
    shell: 
        """
        # Align against reference genome,
        # minimap2 automatically handles
        # conversion of U to T bps 
        minimap2 \\
            -ax splice \\
            -uf \\
            -k14 \\
            {input.ref} \\
            {input.fq} \\
        | samtools view  -Sb - \\
        | samtools sort -@{threads} \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            -
        
        # Gather BAM statistics
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
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    threads: int(allocated("threads", "minimap2_transcriptome", cluster)) 
    shell: 
        """
        # Align against reference genome,
        # minimap2 automatically handles
        # conversion of U to T bps
        minimap2 \\
            -ax splice \\
            -uf \\
            -k14 \\
            {input.ref} \\
            {input.fq} \\
        | samtools view  -Sb - \\
        | samtools sort -@{threads} \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            -
        """
