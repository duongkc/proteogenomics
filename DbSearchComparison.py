#!/usr/bin/python3.6
import sys
import re
import datetime

"""
Short script to compare the PEAKS psm search output to find overlap and distinction between the matched peptides.
"""


def strip_peptide_col(peptide_column):
    """strips peptide column by removing unnecessary information and returns the peptide"""
    sample_line = "K.ELC(+57.02)EQ(+.98)EC(+57.02)EWEEITITGSDGSTR.V"
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep



def parse_genemark_search():
    """Parses the genemark file to extract the peptide column and matching ORF accession"""
    gmfile = "data/propep_g.csv"  # Change to non static later
    with open(gmfile, "r") as f:
        next(f)
        for line in f:
            columns = line.split(",")
            accession = columns[2]
            peptide_col = columns[3]
            peptide = strip_peptide_col(peptide_col)
            print(peptide)


def parse_transdecoder_search():
    """Parses the transdecoder file to extract the peptide column and matching ORFs"""


parse_genemark_search()


def main():
    print("Hello, World!")
