#!/usr/bin/python3
"""
-g data/propep_genemark.csv -t data/propep_transdecoder.csv -p test1
"""


import re

import pandas


def clean_peptide_col(peptide_column):
    """Cleans up the peptide column by removing unnecessary information and returns the peptide"""
    sample_line = "K.ELC(+57.02)EQ(+.98)EC(+57.02)EWEEITITGSDGSTR.V"
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file):
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    csv_data = csv_data[['Protein Accession', 'Peptide']].drop_duplicates(subset=['Peptide'], keep='first')
    print(len(csv_data))
    return csv_data


def find_distinct_peptides(transdecoder_data, genemark_data, prefix):
    """Filters the CSV files so only distinct peptides remain"""
    distinct_td_csv = "comparison_output/{}_distinct_td.csv".format(prefix)
    distinct_gm_csv = "comparison_output/{}_distinct_gm.csv".format(prefix)

    td_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='left', indicator=True) \
        .query("_merge == 'left_only'")
    gm_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='right', indicator=True) \
        .query("_merge == 'right_only'")

    overlap_left = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='left', indicator=True)\
        .query("_merge == 'both'")
    # overlap_merge = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='inner', indicator=True)\
    #     .query("_merge == 'both'")
    # overlap_left[['Protein Accession_x', 'Peptide']]\
    #     .drop_duplicates(subset=['Protein Accession_x', 'Peptide'], keep='first')\
    #     .to_csv(overlap_file, sep=',', mode='w', index=False)
    # print(len(overlap_left.index))
    # print(len(overlap_merge.index))
    print(len(td_merged.index))
    print(len(gm_merged.index))


genemark_file = "data/propep_genemark.csv"
trans_file = "data/propep_transdecoder.csv"

trans_data = extract_csv_data(trans_file)
genemark_data = extract_csv_data(genemark_file)

find_distinct_peptides(trans_data, genemark_data, "test")
