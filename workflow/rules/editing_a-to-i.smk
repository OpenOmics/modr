# A-to-I RNA-editing rules 

# @DISCLAIMER: PLEASE READ THROUGH THIS MESSAGE!
# Rules for running DInoPORE: Direct detection 
# of INOsines in native RNA with nanoPORE sequencing.
# We need to process the data using the exact same set 
# of methods as the authors of DInopore. This will ensure
# the highest level of accuracy and precision when using
# their pre-built models. The DInopore repo on Github is 
# not a general purpose tool. It is hard-coded to run on
# CodeOcean. The paths referenced in the scripts are hard-
# coded to run in the container's filesysem, and there is 
# very little control over each step. Also, if a step fails
# halfway through, it will just re-run everything from the 
# beginning. This is not desirable, but some of these steps 
# can take more than 12 hours to run. The following set of 
# rules below are an attempt to re-engineer the data pro-
# cessing steps described in the DInopore paper/github repo 
# to allow other users to run this awesome tool.
# If you used `modr` or referenced these rules, please 
# do not forget to cite the authors of DInoPORE.
# @Github: 
#    https://github.com/darelab2014/Dinopore
# @Citation:
#   Nguyen, T.A., Heng, J.W.J., Kaewsapsak, P. et al. 
#   Direct identification of A-to-I editing sites with 
#   nanopore native RNA sequencing. Nat Methods 19, 833–844 
#   (2022). https://doi.org/10.1038/s41592-022-01513-3
rule dinopore_graphmap2:
    """
    Data-processing step for setting up dinopore. This step
    converts any U bps to T and aligns the reads to the genome.
    @Input:
        FastQ file (scatter),
        Genomic FASTA file
    @Output:
        Sorted Graphmap2 Genomic BAM file
    """
    input:
        fq  = join(workpath, "{name}", "fastqs", "{name}.fastq.gz"),
        ref = join(workpath, "refs", ref_genome),
    output:
        fq = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq"),
        bam = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sorted.bam"),
        bai = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sorted.bam.bai"),
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
        
        # Convert to bam, index, and sort
        samtools sort -@{threads} \\
            -T "${{tmp}}" \\
            -O bam \\
            -o {output.bam} \\
            {params.sam}
        samtools index -@{threads} \\
            {output.bam} \\
            {output.bai}
        """


rule dinopore_nanopolish:
    """
    Data-processing step align nanopore event signals to k-mers 
    of the reference genome.
    @Input:
        U to T converted FastQ file (scatter),
        Sorted Genomic Graphmap2 BAM file (scatter),
        Genomic FASTA file
    @Output:
        Nanopolish indices,
        Nanopolish signal event align file
    """
    input:
        fq = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq"),
        bam = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sorted.bam"),
        ref = join(workpath, "refs", ref_genome),
    output:
        idx = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq.index"),
        gzi = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq.index.gzi"),
        fai = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq.index.fai"),
        rdb = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.fastq.index.readdb"),
        summary = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.nanopolish.sum"),
        events  = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.nanopolish.eventAlignOut.txt"),
    params:
        rname  = 'dinopolish',
        f5 = lambda w: str(config['fast5'][w.name]),
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_nanopolish", cluster))
    shell: 
        """
        # Build an index mapping
        # basecalled reads to event 
        # signals from the sequencer
        nanopolish index \\
            -d {params.f5} \\
            {input.fq}
        # Overlay event signal to 
        # kmers of reference genome
        nanopolish eventalign \\
            --reads {input.fq} \\
            --bam {input.bam} \\
            --genome {input.ref} \\
            --summary {output.summary} \\
            --print-read-names \\
            --threads {threads} \\
            --scale-events \\
        > {output.events}
        """


rule dinopore_fadict:
    """
    Creates a sequence dictionary of the reference genome 
    and create an index of the genomic FASTA file, which are
    both needed by "sam2tsv.jar".  
    @Input:
        Genomic FASTA file
    @Output:
        Sequence dictionary,
        FASTA index
    """
    input:
        ref = join(workpath, "refs", ref_genome),
    output:
        dct = join(workpath, "refs", "{0}.dict".format(ref_genome)),
        fai = join(workpath, "refs", "{0}.fai".format(ref_genome)),
    params:
        rname  = 'dinodict',
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_fadict", cluster))
    shell: 
        """
        # Build an sequence dictionary
        java -jar ${{PICARDJARPATH}}/picard.jar \\
            CreateSequenceDictionary \\
            R={input.ref} \\
            O={output.dct}
        # Create FASTA index
        samtools faidx \\
            {input.ref}
        """


rule dinopore_sam2tsv:
    """
    Data-processing step to convert bam file into tsv file, then process 
    it to remove (S, H and N), only keeping M (match or mismatch), D (deletion)
    and I (insertion).
    @Input:
        Sorted Graphmap2 Genomic BAM file
        Genomic FASTA file
        Sequence Dictionary
    @Output:
        Filtered Alignments TSV
    """
    input:
        bam = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sorted.bam"),
        ref = join(workpath, "refs", ref_genome),
        dct = join(workpath, "refs", "{0}.dict".format(ref_genome)),
        fai = join(workpath, "refs", "{0}.fai".format(ref_genome)),
    output:
        tsv = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.raw_features.tsv"),
    params:
        rname  = 'dinosam2tsv',
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_sam2tsv", cluster))
    shell: 
        """
        # Convert SAM to raw features TSV
        samtools view \\
            -@{threads} \\
            -h \\
            -F 4 \\
            {input.bam} \\
        | java -jar ${{PICARDJARPATH}}/sam2tsv.jar -r {input.ref} \\
        | awk 'BEGIN{{FS=OFS="\\t"}} ($9 != "S") && ($9 != "H") && ($9 != "N")' - \\
        | awk 'BEGIN{{FS=OFS="\\t"}} \\
            ($7=="."){{$7="-99";}} \\
            ($4=="."){{$4="-99"}} \\
            ($5=="."){{$5="na"}} \\
            ($8=="."){{$8="na"}} \\
            ($9=="D"){{$6=" "}} \\
            ($2==16){{$2="n"}} \\
            ($2==0){{$2="p"}} 1' \\
        > {output.tsv}
        """


rule dinopore_combine_events:
    """
    Data-processing step to combine/collapse event-level signal
    information from nanopolish. For more information, please see:
    https://github.com/darelab2014/Dinopore/blob/main/code/s2.Combine_raw_nnpl.sh
    NOTE: Look into futher optimizing this rule later, some of these 
    calculation can probably be made massively in parallel instead
    of serially.
    @Input:
        Nanopolish signal event align file
    @Output:
        Combined Nanopolish signal event align file
    """
    input:
        events  = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.nanopolish.eventAlignOut.txt"),
    output:
        events  = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.nanopolish.eventAlignOut.combined.txt"),
    params:
        rname  = 'dinocombine',
        outdir = join(workpath, "{name}", "rna-editing", "dinopore"),
        header = join(workpath, "{name}", "rna-editing", "dinopore", 'raw_nanopolish.header'),
        tmp = join(workpath, "{name}", "rna-editing", "dinopore", 'tmp.nnpl'),
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_combine_events", cluster))
    shell: 
        """
        # Setup for the combine step,
        # see github repo for more info:
        # https://github.com/darelab2014/Dinopore/blob/main/code/s2.Combine_raw_nnpl.sh
        cd "{params.outdir}"
        head -n 1 {input.events} \\
            > {params.header}
        noline=$(wc -l {input.events} | cut -d " " -f 1)
        tmp="{params.tmp}"
        num=2
        chunk=10000000
        num1=$(expr $num + $chunk)
        num2=$(expr $num1 + 1)
        filecount=1

        # Combine event-level signal information
        while [ $num -le $noline ]; do
            sed -n "$num,${{num1}}p; ${{num2}}q" {input.events} > "$tmp"
            chk=$(wc -l "$tmp" | cut -f 1 -d " ")
            echo "Processing reads $num to $num1 in file $filecount"
            if [ $chk -gt 0 ]; then
                Rscript ${{DINOPORE_CODE}}/s2.Combine_raw_nanopolish.R \\
                    -f "$tmp" \\
                    -t {threads} \\
                    -o {output.events}.part${{filecount}} \\
                    -s ${{noline}} \\
                    -n ${{num}} \\
                    -c ${{chunk}}
                if test -f tmp.eli; then
                    eli=$(cat tmp.eli)
                    rm "tmp.eli"
                else
                    eli=0
                fi
            fi
            num=$(( $num1 - $eli + 1 ))
            num1=$(expr $num + $chunk)
            num2=$(expr $num1 + 1)
            filecount=$(expr $filecount + 1)
        done

        # Combine results of all parts
        cat ${{DINOPORE_CODE}}/misc/nnpl.header > {output.events}
        cat {output.events}.part* | grep -v contig >> {output.events}
        rm -f tmp.nnpl {output.events}.part*
        """


rule dinopore_feature_table:
    """
    Data-processing step generate a raw feature table for dinopore's 
    pre-built model. For more information, please see the following:
    https://github.com/darelab2014/Dinopore/blob/main/code/S3.Generate_raw_features.sh
    NOTE: Look into futher optimizing this rule later, some of these 
    calculation can probably be made massively in parallel instead
    of serially.
    @Input:
        Filtered Alignments TSV (scatter),
        Combined Nanopolish signal event align file (scatter)
    @Output:
        Combined Nanopolish signal event align file
    """
    input:
        tsv = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.raw_features.tsv"),
        events  = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.nanopolish.eventAlignOut.combined.txt"),
    output: 
        table = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sam2tsv_nanopolish.grp_{group}.txt"),
    params:
        rname  = 'dinofeats',
        outdir = join(workpath, "{name}", "rna-editing", "dinopore"),
        memory = allocated("mem", "dinopore_feature_table", cluster).upper(),
        # List of intermediate files 
        tmp = join(workpath, "{name}", "rna-editing", "dinopore", 'reads.tmp'),
        outname = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.sam2tsv_nanopolish"),
        tsv = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.reads.sam2tsv"),
        events = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.reads.nanopolish"),
        common = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.reads.intersect"),
        intsv  = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.input.filteredSHN.tsv.txt"), 
        innnpl = join(workpath, "{name}", "rna-editing", "dinopore", "{name}.input.nanopolish.eventAlignOut.combined.txt"), 
        # Needs dedicated space for sorting 
        # large files, otherwise /tmp disk
        # quota will get exceeded
        tmpdir = join(workpath, "{name}", "rna-editing", "dinopore", "sort_tmp"), 
    conda: depending(join(workpath, config['conda']['dinopore']), use_conda)
    container: depending(config['images']['dinopore'], use_singularity)
    threads: int(allocated("threads", "dinopore_feature_table", cluster))
    shell: 
        """
        # Setups temporary directory for
        # intermediate files with built-in 
        # mechanism for deletion on exit
        if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
        tmp=$(mktemp -d -p "{params.tmpdir}")
        trap 'rm -rf "${{tmp}}"' EXIT

        # Setup for the feature table step,
        # see github repo for more info:
        # https://github.com/darelab2014/Dinopore/blob/main/code/S3.Generate_raw_features.sh
        cd "{params.outdir}"
        num=1
        chunk=100000
        num1=$(expr $num + $chunk - 1)
        num2=$(expr $num1 + 1)
        filecount=1

        # Extract common reads between filtered
        # alignment TSV and nanopolish events 
        awk -F '\\t' 'NR!=1 {{print $1}}' {input.tsv} \\
            | LC_ALL=C uniq \\
            | LC_ALL=C sort -S {params.memory} --parallel={threads}  \\
        > {params.tsv}
        awk '{{print $2}}' {input.events} \\
            | LC_ALL=C uniq \\
            | LC_ALL=C sort -S {params.memory} --parallel={threads}  \\
        > {params.events}
        # Prints lines/reads common in both files
        comm -12 \\
            {params.tsv} \\
            {params.events} \\
        > {params.common}
        noline=$(wc -l {params.common} | cut -d " " -f 1)

        # Generate raw features files
        while [ $num -le $noline ]; do
            echo "Processing reads $num to $num1 in file $filecount"
            sed -n "$num,${{num1}}p; ${{num2}}q" {params.common} \\
            > {params.tmp}
            chk=$(wc -l {params.tmp} | cut -f 1 -d " ")
            if [ $chk -gt 0 ]; then
                # The headers can also be found here:
                # https://github.com/darelab2014/Dinopore/tree/main/code/misc
                cat ${{DINOPORE_CODE}}/misc/nnpl.header \\
                > {params.innnpl}
                cat ${{DINOPORE_CODE}}/misc/tsv.header \\
                > {params.intsv}
                LC_ALL=C grep \\
                    -h \\
                    -F \\
                    -w \\
                    -f {params.tmp} \\
                    {input.tsv} \\
                >> {params.intsv}
                LC_ALL=C grep \\
                    -h \\
                    -F \\
                    -w \\
                    -f {params.tmp} \\
                    {input.events} \\
                >> {params.innnpl}

                # Create chunked raw features files,
                # Rscript cannot use files with abs paths
                rel_intsv="$(basename {params.intsv})"
                rel_innnpl="$(basename {params.innnpl})"
                Rscript ${{DINOPORE_CODE}}/s3.Generate_raw_feature_table.R \\
                    -n ${{rel_innnpl}} \\
                    -t ${{rel_intsv}} \\
                    -o {params.outname}.part${{filecount}}

                num=$num2
                num1=$(expr $num + $chunk - 1)
                num2=$(expr $num1 + 1)
                filecount=$(expr $filecount + 1)
                rm {params.tmp}
            fi
        done

        # Combine results of all parts
        cat ${{DINOPORE_CODE}}/misc/inae.header \\
        > {output.table}
        cat {params.outname}.part*.positive \\
            | LC_ALL=C sort -T "${{tmp}}" -k1,1 -k3,3n \\
            | grep -v contig \\
        >> {output.table}
        cat {params.outname}.part*.negative \\
            | LC_ALL=C sort -T "${{tmp}}" -k1,1 -k3,3n \\
            | grep -v contig \\
        >> {output.table}
        
        # Cleanup intermediate files
        rm -f {params.tsv} \\
            {params.events} \\
            {params.common} \\
            {params.intsv} \\
            {params.innnpl} \\
            {params.outname}.part*.positive \\
            {params.outname}.part*.negative
        """