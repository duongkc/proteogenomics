#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.
"""


import re
import sys

import pandas
from Bio import SeqIO


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
    csv_data = csv_data.drop_duplicates(subset='Peptide', keep='first').reset_index(drop=True)
    return csv_data


def search_peptide_db(peptide_data, database):
    """Checks for presence of peptides in protein database"""
    # orig = 0
    # counter = 0
    for i, row in peptide_data.iterrows():
        print(i)
    #     orig += 1
    #     flag = 0
    #     peptide = row['Peptide']
    #     for record in SeqIO.parse(database, "fasta"):
    #         if peptide in record.seq:
    #             flag = 1
    #     if not flag:
    #         counter += 1
    # print(counter)
    # print(orig)


def main():
    csv_data = extract_csv_data("data/propep_g.csv")
    database_file = "data/sample_sprot.fasta"
    with open(database_file, "r") as database:
        search_peptide_db(csv_data, database)


if __name__ == '__main__':
    main()
