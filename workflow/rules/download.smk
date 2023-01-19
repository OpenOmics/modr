# Download reference rules
rule download_genome:
    """
    Setup step to download references files for the pipeline.
    This step downloads the genomic fasta file from GENOCODE.  
    @Input:
        None
    @Output:
        Genomic fasta file
    """
    output:
        genome = join(workpath, "refs", "{ref}.genome.fa"),
    params:
        rname = 'getgenome',
        # Download links to each of the 
        # ref files from GENCODES FTP
        uri = config['references'][genome]['GENOME_FA'],
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    shell: 
        """
        # Download primary assembly,
        # reference genome sequence
        wget {params.uri}
        gzip -d {output.genome}.gz
        rm {output.genome}.gz
        """


rule download_transcriptome:
    """
    Setup step to download references files for the pipeline.
    This step downloads the transcriptomic fasta file from GENOCODE.  
    @Input:
        None
    @Output:
        Transcriptomic fasta file
    """
    output:
        transcripts = join(workpath, "refs", "{ref}.transcripts.fa"),
    params:
        rname = 'gettranscipts',
        # Download links to each of the 
        # ref files from GENCODES FTP
        uri = config['references'][genome]['TRANSCRIPTS_FA'],
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    shell: 
        """
        # Download transcript sequences
        wget {params.uri}
        gzip -d {output.transcripts}.gz
        rm {output.transcripts}.gz
        """


rule download_gtf:
    """
    Setup step to download references files for the pipeline.
    This step downloads the annotation file (GTF) from GENOCODE.  
    @Input:
        None
    @Output:
        GTF file
    """
    output:
        gtf = join(workpath, "refs", "{ref}.gtf"),
    params:
        rname = 'getgtf',
        # Download links to each of the 
        # ref files from GENCODES FTP
        uri = config['references'][genome]['GTF_FILE'],
    conda: depending(join(workpath, config['conda']['modr']), use_conda)
    container: depending(config['images']['modr'], use_singularity)
    shell: 
        """
        # Download PRI gene annotation
        wget {params.uri}
        gzip -d {output.gtf}.gz
        rm {output.gtf}.gz
        """
