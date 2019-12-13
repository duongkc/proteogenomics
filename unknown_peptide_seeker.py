#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.
"""


import re
import sys

import pandas


def clean_peptide_col(peptide_column):
    """Cleans up the peptide column by removing unnecessary information and returns the peptide"""
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file):
    """Reads PEAKS protein-peptide.csv file as dataframe to extract unique peptides and matching ORF accessions"""
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    csv_data = csv_data.drop_duplicates(subset=['Protein Accession', 'Peptide'], keep='first')
    return csv_data