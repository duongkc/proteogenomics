#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.
"""
import datetime
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

    flag_list = [0] * len(peptide_data.index)
    for record in SeqIO.parse(database, "fasta"):
        for i, row in peptide_data.iterrows():
            peptide = row['Peptide']
            if peptide in record.seq:
                flag_list[i] = 1
    print(flag_list)


def main():
    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    csv_data = extract_csv_data("data/propep_g.csv")
    database_file = "data/uniprot_sprot.fasta"
    with open(database_file, "r") as database:
        search_peptide_db(csv_data, database)
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main()
