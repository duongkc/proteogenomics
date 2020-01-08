#!/usr/bin/python3
"""
A short module that converts PEAKS protein-peptide.csv files to pandas dataframes in which the peptide column has
been tidied up for further analysis purposes.
"""
import re

import pandas


def clean_peptide_col(peptide):
    """Cleans up the peptide column entry by removing unnecessary information and returns the peptide"""
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file):
    """Reads PEAKS protein-peptide.csv file as dataframe with only unique peptides present"""
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    csv_data = csv_data.drop_duplicates(subset=['Peptide'], keep='first').reset_index(drop=True)
    return csv_data
