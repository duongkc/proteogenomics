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


def parse_csv(input_file, db_type):
    """Parses the PEAKS protein-peptide csv file to extract the peptides and matching ORF accessions"""
    # And write them to a new temporary file, (should be split into sub methods tbh
    if db_type == "genemark":
        temp_output = "output/temp_gm_output.csv"
        temp_just_peptides = "output/temp_gm_just_peps.csv"
    else:
        temp_output = "output/temp_td_output.csv"
        temp_just_peptides = "output/temp_td_just_peps.csv"

    with open(temp_output, "w+") as new_file, open(temp_just_peptides, "w+") as temp_file:
        new_file.write("Protein Accession, Peptide\n")
        with open(input_file, "r") as f:
            next(f)
            for line in f:
                columns = line.split(",")
                accession = columns[2]
                peptide_col = columns[3]
                peptide = strip_peptide_col(peptide_col)
                new_line = accession + "," + peptide + "\n"
                new_file.write(new_line)
                temp_file.write(peptide + "\n")


def parse_transdecoder_search():
    """Parses the transdecoder file to extract the peptide column and matching ORFs"""


parse_csv("data/propep_g.csv", "genemark")
parse_csv("data/propep_t.csv", "transdecoder")


def main():
    print("Hello, World!")
