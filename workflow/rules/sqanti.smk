# Sqanti related quality-control and filtering rules,
# Sqanti is being used to annotate/characterize novel 
# isoforms and to build an even higher-confidence,
# filtered set of unique transcripts from flair.
# The resulting annotation/transcriptome will be
# used to quantify known/novel isoforms.
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
    pass


rule sqanti_ml_filter:
    """
    Data-processing step to filter the sqanti qc output. The auhtor
    from sqanti highly recommends filtering its output before using
    it in down-stream analysis. Sqanti has a new filtering method
    that employs random forest to discriminate potential artifacts
    from true isoforms without the need for user-defined rules or
    manually-set thresholds (i.e. previous method). For more info, 
    please read through sqanti3's documenation:
    https://github.com/ConesaLab/SQANTI3/wiki/
    Github: https://github.com/ConesaLab/SQANTI3
    @Input:
        Sqanti Classification file (TSV) 
    @Output:
        ML Filtered Sqanti Classification file (TSV),
        ML Filtered Corrected Annotation (GTF),
        ML Filtered Corrected Transcriptome (FASTA)
    """
    pass