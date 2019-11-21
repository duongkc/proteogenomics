#!/usr/bin/python3.6
import sys
import datetime

"""
Short script to compare the PEAKS psm search output to find overlap and distinction between the matched peptides.
"""


def extract_peptide_column():
    """extracts peptide column and removes unnecessary information, such as the modification between parentheses"""


def parse_genemark_search():
    """Parses the genemark file to extract the peptide column and matching ORF"""
    gmfile = "data/genemark_dbsearch_psm.csv"  # Change to non static later
    with open(gmfile, "r") as f:
        next(f)
        for line in f:
            columns = line.split(",")
            print(columns[0])


def parse_transdecoder_search():
    """Parses the transdecoder file to extract the peptide column and matching ORFs"""


parse_genemark_search()


def main():
    print("Hello, World!")
