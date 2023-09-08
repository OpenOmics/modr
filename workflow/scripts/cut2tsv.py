#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
import sys, os
import textwrap

# 3rd party packages,
# installed from pypi
import argparse
import pandas as pd


_help = textwrap.dedent("""
@Usage:
    $ ./cut2tsv.py [-h] \\
         -c COL_A [COL_B COL_Y COL_Z] \\
         -i FILE

@About:
    Given an input file and 1 or more column
    names to extract, this script will parse
    each of those columns by their name. The 
    script can also be used to re-order the
    columns of a file by their name too. With
    that being said, all operations performed
    by this script work on column names.
    
    To perform specific operations on columns
    by their column their index (i.e. 1, 2, 3
    ...), please use the built-in cut command
    or awk.

    All output is directed to standard output
    in (awk-friendly) tab-delimited format. It
    can be captured, redirected, or piped into
    another process. 

@Required Arguments:
    -c, --column COL  One or more column names
                      to extract. The order of
                      the columns will be the 
                      the order of the output.
                      If a column is provided
                      that does not exist in 
                      the input file, then a 
                      warning message will be
                      produced, and it will 
                      skip over the column. 

    -i, --input FILE  Input file to extract or
                      reorder column(s) by their
                      name. The following input 
                      file types are types are
                      supported:
                        + TSV (tab-delimited)
                        + CSV (comma-delimited)
                        + Excel
@Options
    --h, --help     Shows help message and exits.


@Example
    # Input file
    $ cat file.tsv
    a   b   c
    1   2   3
    4   5   6
    7   8   9

    # Subset and re-order
    # the example input file 
    $ cut2tsv.py -c b a -i file.tsv
    b   a
    2   1
    5   4
    8   7
"""
)


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


def read(filename, subset=[], skip='#', **kwargs):
    """Reads in an input file as a dataframe. Determines the 
    correct handler for reading in a given MAF file. Supports reading
    in TSV files (.tsv, .txt, .text, .vcf, or .maf), CSV files (.csv), 
    and excel files (.xls, .xlsx, .xlsm, .xlsb, .odf, .ods, .odt ). 
    The subset option allows a users to only select a few columns 
    given a list of column names.
    @param filename <str>:
        Path of an MAF-like file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    # Get file extension
    extension = os.path.splitext(filename)[-1].lower()

    # Assign a handler to read in the file
    if extension in ['.xls', '.xlsx', '.xlsm', '.xlsb', '.odf', '.ods', '.odt']:
        # Read in as an excel file
        return excel(filename, subset, skip, **kwargs)
    elif extension in ['.csv']:
        # Read in as an CSV file
        return csv(filename, subset, skip, **kwargs)
    else:
        # Default to reading in as an TSV file
        # Tab is the normal delimeter for MAF or VCF files
        # MAF files usually have one of the following
        # extensions: '.tsv', '.txt', '.text', '.vcf', '.maf'
        return tsv(filename, subset, skip, **kwargs)


def excel(filename, subset=[], skip='#', **kwargs):
    """Reads in an excel file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
        Path of an EXCEL file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset: 
        return pd.read_excel(filename, comment=skip, **kwargs)[subset]
    
    return pd.read_excel(filename, comment=skip, **kwargs)


def tsv(filename, subset=[], skip='#', **kwargs):
    """Reads in an TSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
        Path of an TSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset: 
        return pd.read_table(filename, comment=skip, **kwargs)[subset]
    
    return pd.read_table(filename, comment=skip, **kwargs)


def csv(filename, subset=[], skip='#', **kwargs):
    """Reads in an CSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of 
    column names.
    @param filename <str>:
        Path of an CSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset: 
        return pd.read_csv(filename, comment=skip, **kwargs)[subset]
    
    return pd.read_csv(filename, comment=skip, **kwargs)

def parse_arguments():
    """Returns an argparse object of parsed command-line args."""
	# Parse command-line arguments
    # Create a top-level parser
    parser = argparse.ArgumentParser(
        usage = argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = _help,
        add_help=False
    )

    # Required Positional Arguments
    # List of input files
    parser.add_argument(
        '-c', '--column',
        required=True,
        nargs = '+',
        help = argparse.SUPPRESS
    )
    # Output file name
    parser.add_argument(
        '-i', '--input',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )

    # Options
    # Add custom help message
    parser.add_argument(
        '-h', '--help', 
        action='help', 
        help=argparse.SUPPRESS
    )

    # Collect parsed arguments
    args = parser.parse_args()

    return args


def main():
    """Collect command line args and build the matrix."""
	# Parse command-line arguments
    args = parse_arguments()

    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: cut2tsv.py [-h] ...')

    # Read in first line of
    # input file to get its
    # column names 
    input_header = read(
        args.input,
        nrows=1
    ).columns.tolist()

    # Check if user columns do not
    # exist in the input file
    extract_columns = args.column
    header_set = set(input_header)
    extract_set = set(extract_columns)
    missing = extract_set - header_set
    if missing:
        # Display warning to stderr 
        # and continue with processing 
        err(
            "Warning: The following columns do not " +
            "exist in your input file: '{}'. ".format(missing) +
            "Skipping over these columns..."
        )
        # Remove missing columns
        extract_columns = [c for c in extract_columns if c not in missing]

    # Check if there is any to extract
    # or re-order (i.e. non-empty list)
    if not extract_columns:
        fatal(
            "Error: None of the provided column names " +
            "exist in your input file!"
        )

    # Subset and/or re-order columns
    df = read(
        args.input,
        subset = extract_columns
    )

    # Write the resulting df to stdout
    df.to_csv(
        sys.stdout, 
        sep = "\t", 
        header = True, 
        index = False
    )


if __name__ == '__main__':
    # Create the matrix
    try:
        main()
    except BrokenPipeError:
        pass