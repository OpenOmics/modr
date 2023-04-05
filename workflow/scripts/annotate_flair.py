#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
ABOUT: Annotate FLAIR's results with gene/transcript information.
REQUIRES:
  - python>=3
DISCLAIMER:
                    PUBLIC DOMAIN NOTICE
        NIAID Collaborative Bioinformatics Resource (NCBR)
   National Institute of Allergy and Infectious Diseases (NIAID)
This software/database is a "United  States Government Work" under
the terms of the United  States Copyright Act.  It was written as 
part of the author's official duties as a United States Government
employee and thus cannot be copyrighted. This software is freely
available to the public for use.
Although all  reasonable  efforts have been taken  to ensure  the
accuracy and reliability of the software and data, NCBR do not and
cannot warrant the performance or results that may  be obtained by 
using this software or data. NCBR and NIH disclaim all warranties,
express  or  implied,  including   warranties   of   performance, 
merchantability or fitness for any particular purpose.
Please cite the author and NIH resources like the "Biowulf Cluster" 
in any work or product based on this material.
"""

# Python standard library
from __future__ import print_function
import sys

# Additional script metadata
__version__ = 'v0.1.0'
__authors__ = 'Skyler Kuhn'
__email__   = 'skyler.kuhn<AT>nih.gov'

usage = \
"""
USAGE:
  $ ./annotate.py [-h] \\
      <flair.transcripts.tpm.tsv> \\
      <gencode.v42.primary_assembly.annotation.gtf> 

SYNOPSIS:
  Annotates each recorded FLAIR transcript. Each transcript
  will be anntotated as either known or novel. In addition,
  the following information will be added: 
    - 'location'
    - 'gene_name'
    - 'transcript_name'
    - 'biotype' 
    - 'ensembl_url'

   Please note that the order of the positional command-line
   arguments matter. The first positional argument provided 
   to this script should be the FLAIR's counts (raw or tpm),
   while the second positional argument is the GTF file that
   was used for correction (i.e. `flair correct --gtf` file). 

   The annotated counts are returned to standard output.

EXAMPLE:
  $ ./annotate.py \\
      flair.transcripts.count.tsv \\
      gencode.v42.primary_assembly.annotation.gtf \\
    > flair.transcripts.count.annotated.tsv

VERSION:
  {0}
""".format(__version__)


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def clean(string):
    """Cleans a string to remove double and single quotes."""
    return string.replace('"', '').replace("'", '')


def extract_9th_column_information(key_value_line, values_to_extract):
    """Given a string containing key/value information in
    the 9th column in the GTF file, this will extract info
    for a list of known values."""
    kv_list = key_value_line.lstrip().rstrip().rstrip(';').split(';')
    metadata = {attribute: '' for attribute in values_to_extract}
    for kv in kv_list:
        k,v = kv.lstrip().split(' ', maxsplit=1)
        k = clean(k)
        if k in values_to_extract:
            metadata[k] = v
    return metadata


def default_lookup(lookup):
    """Sets a default value for a dictionary look."""
    try:
        v = lookup
    except KeyError:
        v = ''
    return v


def parsed_gtf(gencode_gtf_file):
    """Parses a GENCODE GTF file to get additional metadata
    for each annotated transcript. YMMV for GTF files from 
    providers outside of GENCODE.
    @param gencode_gtf_file <str>: 
        GENCODE GTF file to parse
    @returns gtf_metadata <dict[ENSG|ENST]={key_9th_col:value}>
    """
    gtf_metadata = {}

    with open(gencode_gtf_file) as gtf:
        for line in gtf:
            linelist = line.lstrip().rstrip().split('\t')
            additional_information = {}
            if line.startswith('#'):
                continue    # skip over comment lines

            # Parse important information for all the 
            # annotated genes and transcripts
            feature = linelist[2].lower()
            if feature == 'gene' or feature == 'transcript':
                # Parse important information from the 9th 
                # column which contains key,value pairs of 
                # metadata for a given gene/transcript 
                chr_start_stop = "{0}:{1}-{2}".format(
                    linelist[0],
                    linelist[3],
                    linelist[4]
                )
                if feature == 'gene':
                    # Genes contains a subset of 
                    # features we want  
                    nineth_column_information = extract_9th_column_information(
                        linelist[8], 
                        [
                            'gene_id',
                            'gene_type',
                            'gene_name'
                        ] 
                    )
                    # Create lookup for ENSG* information
                    feature_id = clean(nineth_column_information['gene_id'])

                if feature == 'transcript':
                    nineth_column_information = extract_9th_column_information(
                        linelist[8],
                        [
                            'gene_id',
                            'transcript_id',
                            'gene_type',
                            'gene_name',
                            'transcript_name',
                            'transcript_type'
                        ] 
                    )
                    # Create lookup for ENST* information
                    feature_id = clean(nineth_column_information['transcript_id'])

                # Add feature's annotation information                 
                gtf_metadata[feature_id] = nineth_column_information
                gtf_metadata[feature_id]['chr_start_stop'] = chr_start_stop
                # Add some default annotation information,
                # for when lookups fail OR we come across
                # an edge-case not currently supported.
                gtf_metadata['__default__'] = {
                    'chr_start_stop': '',
                    'gene_id': '',
                    'transcript_id': '',
                    'gene_type': '',
                    'gene_name': '',
                    'transcript_name': '',
                    'transcript_type': ''
                }
    
    return gtf_metadata


def main():
    # Check for help option and exit
    if '-h' in sys.argv or '--help' in sys.argv:
        err(usage)
        sys.exit()
    try:
        # Parses command line arguments,
        # and checks for basic usage.
        # Please NOTE: the order of 
        # positional args matter!
        flair_results = sys.argv[1]
        gtf_file = sys.argv[2]
    except IndexError: 
        # Missing required args
        err('Error: missing one or more required argument(s)!')
        fatal(usage)

    # Get annotation information from
    # the GENCODE GTF file
    gtf_annotations = parsed_gtf(gtf_file)
    # Add parsed annotation information to results,
    # FLAIR reports known and novel isoforms as follows:
    # 1. Known ~ ENST00000426218.1_ENSG00000229271.3
    #      The id of these transcripts are just a composite
    #      of the Ensembl transcript_id and gene_id.
    # 2. Novel ~ 17a55928-3b95-4d82-8c1f-9d8187389cd5_ENSG00000241837.7
    #      If the isoform has a splice pattern in common with 
    #      a gene but doesn't match any of the known transcripts, 
    #      a randomly selected read name is used in combination 
    #      with the gene_id.
    # 3. Novel ~ 87ab9469-20d0-44e9-8d61-6182f919aeb9_chr21:33903000
    #      Lastly, if an isoform does not match a gene, a randomly 
    #      selected read name is combined with the genome location. 
    with open(flair_results, 'r') as results:
        header = next(results).lstrip().rstrip().split('\t')
        # Insert new information after
        # the first column (ids)
        header[1:1] = [
            'location',
            'gene_name',
            'transcript_name',
            'biotype',
            'ensembl_url'
        ]; print('\t'.join(header))
        
        for line in results:
            # Depending on whether transcript is known or novel, the 
            # following additional metadata will be added:
            #    - 'chr_start_stop'
            #    - 'gene_name'
            #    - 'transcript_name'
            #    - 'biotype' 
            #    - 'EnsemblLink'
            linelist = line.lstrip().rstrip().split('\t', maxsplit=1)
            flair_id = linelist[0].split('_')

            if flair_id[0].startswith('ENST'):
                # 1. Known transcript: Identifer is composite of
                #    transcript_id + gene_id, see examples below: 
                #       - ENST00000426218.1_ENSG00000229271.3
                #       - ENST00000706247.1-0_ENSG00000185829.19
                feature_id = flair_id[0].split('-')[0]    # second example above
                try:
                    features = gtf_annotations[feature_id]
                except KeyError:
                    # Cannot find information for feature,
                    # set to default values (empty strings)
                    err('Warning cannot find information for "{0}"'.format(feature_id))
                    features = gtf_annotations['__default__']
                location  = features['chr_start_stop']
                gene_name = features['gene_name']
                transcript_name = features['transcript_name']
                biotype = features['transcript_type']
                url = "https://useast.ensembl.org/Gene/Summary?g={0}".format(
                    flair_id[-1].split('.', maxsplit=1)[0]
                )

            elif flair_id[-1].startswith('ENSG'):
                # 2. Novel transcript in gene: Identifer contains
                #    the gene_id, see two examples below: 
                #       - 17a55928-3b95-4d82-8c1f-9d8187389cd5_ENSG00000241837.7
                #       - a7b8140f-cd36-4f40-ae6a-583d73d467b2-0_ENSG00000169084.15-PAR-Y
                feature_id = flair_id[-1].replace('-','_')      # PAR-Y needs to be PAR_Y
                try:
                    features = gtf_annotations[feature_id]
                except KeyError:
                    err('Warning cannot find information for "{0}"'.format(feature_id))
                    features = gtf_annotations['__default___']  # set to empty strings
                location = features['chr_start_stop']
                gene_name = features['gene_name']
                transcript_name = ''
                biotype = gtf_annotations[feature_id]['gene_type']
                url = "https://useast.ensembl.org/Gene/Summary?g={0}".format(
                    flair_id[-1].split('.', maxsplit=1)[0]
                )

            else:
                # 3. Novel transcript, not in known gene: Identifer 
                #    constists of chr:start position, see example: 
                #       -  87ab9469-20d0-44e9-8d61-6182f919aeb9_chr21:33903000
                position = int(flair_id[-1].split(':')[-1])
                start = max(position-500, 1)
                stop = position + 500
                location = "{0}:{1}-{2}".format(
                    flair_id[-1].split(':')[0],
                    start,
                    stop
                ) # chr_start_stop
                gene_name = ''
                transcript_name = ''
                biotype = ''
                url = ''
            
            # Add additional metadata from the GTF file 
            additional_metadata = [
                location,
                gene_name,
                transcript_name,
                biotype,
                url
            ]
            linelist[1:1] = additional_metadata
            print('\t'.join(linelist))


if __name__ == '__main__':
    # Call main method 
    main()