# A-to-I RNA-editing rules 

# Rules for running DInoPORE: Direct detection 
# of INOsines in native RNA with nanoPORE sequencing.
# We need to process the data using the exact same set 
# of methods as the authors of DInopore. This will ensure
# the highest level of accuracy and precision when using
# their pre-built models. The DInopore repo on Github is 
# not a general purpose tool. It is hard-coded to run on
# CodeOcean. The paths referenced in the scripts are hard-
# coded to run in the container's filesysem and there is 
# very little control over each step. Also, if a step fails
# halfway through, it will just re-run everything from the 
# beginning. This is not desirable, but some of these steps 
# can take more than 8 hours to run. The following set of 
# rules below are an attempt to re-engineer the data pro-
# cessing steps described in the DInopore paper/github repo 
# to allow other users to run this awesome tool.  
# @Github: 
#    https://github.com/darelab2014/Dinopore
# @Citation:
#   Nguyen, T.A., Heng, J.W.J., Kaewsapsak, P. et al. 
#   Direct identification of A-to-I editing sites with 
#   nanopore native RNA sequencing. Nat Methods 19, 833–844 
#   (2022). https://doi.org/10.1038/s41592-022-01513-3
# If you used `modr` or referenced these rules, please 
# do not forget to cite the authors of DInoPORE.
rule dinopore_graphmap2:
    """
    Data-processing step for setting up dinopore. This step
    converts any U bps to T and aligns the reads to the genome.
    @Input:
        Nanofilt quality filtered FastQ file (scatter),
        Genomic FASTA file
    @Output:
        Graphmap2 Genomic BAM file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.fastq.gz"),
        ref = join(workpath, "refs", ref_genome),
    output:
        fq = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq"),
        bam = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.bam"),
        bai = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.bam.bai"),
    params:
        rname  = 'dinomap',
        sam = join(workpath, "{name}", "rna-editing", "dinopore", "genome_tmp", "{name}.sam"),
        tmpdir = join(workpath, "{name}", "rna-editing", "dinopore", "genome_tmp"),
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_graphmap2", cluster)) 
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Convert any U bps to T 
        zcat {input.fq} \\
            | sed 's/U/T/g' \\
        >  {output.fq}

        # Align against genome with graphmap2
        graphmap2 align \\
            -x rnaseq \\
            -t {threads} \\
            -r {input.ref} \\
            -d {output.fq} \\
            -o {params.sam}
        
        # Convert to bam, index and sort
        samtools sort -@{threads} \\
            -T "${{tmp}}" \\
            -O bam \\
            --write-index \\
            -o {output.bam}##idx##{output.bai} \\
            {params.sam}
        """
