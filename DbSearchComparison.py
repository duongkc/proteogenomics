#!/usr/bin/python3
import sys
import datetime
import re
import pandas

"""
Short script to compare the PEAKS psm search output to find overlap and distinction between the matched peptides.
"""


def clean_peptide_col(peptide_column):
    """Cleans up the peptide column by removing unnecessary information and returns the peptide"""
    sample_line = "K.ELC(+57.02)EQ(+.98)EC(+57.02)EWEEITITGSDGSTR.V"
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file, db_type):
    """Parses the PEAKS protein-peptide csv file to extract the peptides and matching ORF accessions"""
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    csv_data = csv_data.drop_duplicates(subset=['Protein Accession', 'Peptide'], keep='first')
    return csv_data


def find_distinct_peptides(transdecoder_data, genemark_data):
    """Filters the CSV files so only distinct peptides remain"""
    with open("output/distinct_td.csv", "w+") as distinct_transdecoder, \
            open("output/distinct_gm.csv", "w+") as distinct_genemark:
        td_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='left', indicator=True) \
            .query("_merge == 'left_only'")
        gm_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='right', indicator=True) \
            .query("_merge == 'right_only'")
        td_merged[['Protein Accession_x', 'Peptide']] \
            .to_csv(distinct_transdecoder, sep=',', mode='w', index=False, header=['Protein Accession', 'Peptide'])
        gm_merged[['Protein Accession_y', 'Peptide']] \
            .to_csv(distinct_genemark, sep=',', mode='w', index=False, header=['Protein Accession', 'Peptide'])


genemark_data = extract_csv_data("data/propep_g.csv", "genemark")
transdecoder_data = extract_csv_data("data/propep_t.csv", "transdecoder")

find_distinct_peptides(transdecoder_data, genemark_data)


def main():
    print("Hello, World!")
