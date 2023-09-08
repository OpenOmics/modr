#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
import sys

# Local imports
from utils import (
    Colors,
    err,
    fatal
)


def clean(s, remove=['"', "'"]):
    """Cleans a string to remove any defined leading or trailing characters.
    @param s <str>:
        String to clean.
    @param remove list[<str>]:
        List of characters to remove from beginning or end of string 's'.
    @return s <str>:
        Cleaned string
    """
    for c in remove:
        s = s.strip(c)
    return s


def contrasts(file, groups, delim='\t'):
    """Reads and parses the group comparison file, contrasts.tsv, into a 
    dictionary. This file acts as a config file to setup contrasts between
    two groups, where groups of samples are defined in the groups.tsv file.
    This information is used in differential analysis, like differential 
    gene expression, etc. 
    @Example: contrasts.tsv
        G2  G1
        G4  G3
        G5  G1
    >> contrasts = contrasts('contrasts.tsv', groups = ['G1', 'G2', 'G3', 'G4', 'G5'])
    >> contrasts
    [
        ["G2",  "G1"],
        ["G4",  "G3"],
        ["G5",  "G1"]
    ]
    @param file <str>:
        Path to contrasts TSV file.
    @param groups list[<str>]:
        List of groups defined in the groups file, enforces groups exist.
    @return comparisons <list[list[str, str]]>:
        Nested list contain comparsions of interest.  
    """

    c = Colors()
    errors = []
    comparsions = []
    line_number = 0
    with open(file) as fh:
        for line in fh:
            line_number += 1
            linelist = [clean(l.strip()) for l in line.split(delim)]
            try:
                g1 = linelist[0]
                g2 = linelist[1]
                if not g1 or not g2: continue # skip over empty lines
            except IndexError:
                # Missing a group, need two groups to tango
                # This can happen if the file is NOT a TSV file,
                # and it is seperated by white spaces, :(  
                err(
                '{}{}Warning: {} is missing at least one group on line {}: {}{}'.format(
                    c.bg_yellow,
                    c.black,
                    file,
                    line_number,
                    line.strip(),
                    c.end
                    )
                )
                err('{}{}\t  └── Skipping over line, check if line is tab seperated... {}'.format(
                    c.bg_yellow,
                    c.black,
                    c.end)
                )
                continue
            # Check to see if groups where defined already,
            # avoids user errors and spelling errors
            for g in [g1, g2]:
                if g not in groups:
                    # Collect all error and report them at end
                    errors.append(g)
            
            # Add comparsion to list of comparisons
            if [g1, g2] not in comparsions:
                comparsions.append([g1, g2])

    if errors:    
        # One of the groups is not defined in groups file
        err('{}{}Error: the following group(s) in "{}" are not defined in --groups file! {}'.format(
            c.bg_red, 
            c.white,
            file,
            c.end)
        )
        fatal('{}{}\t  └── {} {}'.format(
            c.bg_red,
            c.white,
            ','.join(errors),
            c.end)
        )
    
    return comparsions


if __name__ == '__main__':
    # Testing TSV parser
    groups = {"T1": ["A","B"], "T2": ["C","D"]}
    comparsions = contrasts(sys.argv[2], groups=groups.keys())
    print(comparsions)