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


def join_dataframes(data):
    joined_dataframe = pd.DataFrame()
    with open(data, "r") as file_list:
        for file in file_list:
            csv_data = extract_csv_data(file.strip(), drop_dupes=True)
            joined_dataframe = joined_dataframe.append(csv_data[['Peptide']], ignore_index=True)
    return joined_dataframe


def create_peptide_list(left_file, right_file):
    """Creates a list of all peptides as a DataFrame column"""
    joined_left = join_dataframes(left_file)
    joined_right = join_dataframes(right_file)
    all_peptides = joined_left.append(joined_right, ignore_index=True) \
        .drop_duplicates(subset=['Peptide'], keep='first').reset_index(drop=True)
    with open("output/all_peptides_unknown_td.csv", "w+") as output:
        all_peptides.to_csv(output, sep=',', mode='w', line_terminator='\n')


def trim_first_last(peptide_file):
    """Removes first and last character from peptide sequence"""
    data = pd.read_csv(peptide_file, header='infer', delimiter=',', index_col=0)
    for i, row in data.iterrows():
        row['Peptide'] = row['Peptide'][1:-1]
    with open("output/all_peptides_unknown_gm_trim.csv", "w+") as output:
        data.to_csv(output, sep=',', mode='w', line_terminator='\n')
