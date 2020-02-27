#!/usr/bin/python3
"""
A short module that converts PEAKS protein-peptide.csv files to pandas dataframes in which the peptide column has
been tidied up for further analysis purposes.
"""
import re

import pandas as pd


def clean_peptide_col(peptide):
    """Cleans up the peptide column entry by removing unnecessary information and returns the peptide"""
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file, drop_dupes):
    """Reads PEAKS protein-peptide.csv file as dataframe with only unique peptides present"""
    csv_data = pd.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    if drop_dupes:
        csv_data = csv_data.drop_duplicates(subset=['Peptide'], keep='first').reset_index(drop=True)
    return csv_data


def trim_first_last(peptide_file):
    """Removes first and last character from peptide sequence"""
    data = pd.read_csv(peptide_file, header='infer', delimiter=',', index_col=0)
    for i, row in data.iterrows():
        row['Peptide'] = row['Peptide'][1:-1]
    with open("output/all_peptides_unknown_gm_trim.csv", "w+") as output:
        data.to_csv(output, sep=',', mode='w', index=False, line_terminator='\n')
