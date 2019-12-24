#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.
"""
import datetime
import gzip
import multiprocessing as mp
import re
import sys

import pandas
from Bio import SeqIO
import numpy as np


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


def flag_peptides(record, peptide_data):
    peptide_data.reset_index(drop=True)
    flag_list = [0] * len(peptide_data.index)
    for i, row in peptide_data.iterrows():
        peptide = row['Peptide']
        if peptide in record.seq:
            flag_list[i] = 1
    return flag_list


def search_peptide_db(arguments):
    """Checks for presence of peptides in protein database"""
    n = 4
    peptide_data, database_file, offset = arguments

    with open(database_file, "r") as database:
        count = 1
        flag_list = [0] * len(peptide_data.index)  # REVERSE THIS
        for record in SeqIO.parse(database, "fasta"):
            if (count - offset) % n == 0:
                for i, row in peptide_data.iterrows():
                    peptide = row['Peptide']
                    if peptide in record.seq:
                        flag_list[i] = 1
            count += 1
    return flag_list


def write_unknown_peptide_data(peptide_data, flags):
    output = "output/test_unknowns.csv"
    flag_list_0 = np.invert(np.array(flags[0], dtype=bool))
    flag_list_1 = np.invert(np.array(flags[1], dtype=bool))
    flag_list_2 = np.invert(np.array(flags[2], dtype=bool))
    flag_list_3 = np.invert(np.array(flags[3], dtype=bool))

    # print(flag_list_0)
    # print(flag_list_1)
    # print(flag_list_2)
    # print(flag_list_3)
    # print(peptide_data.iloc[3])
    # print(peptide_data.iloc[6])
    # print(peptide_data.iloc[24])
    # print(peptide_data.iloc[48])
    # print(peptide_data.iloc[58])
    # print(peptide_data.iloc[123])
    merged_flag_list = flag_list_0 & flag_list_1 & flag_list_2 & flag_list_3
    # print(merged_flag_list)
    with open(output, "w") as unknown_pep_file:
        filtered_df = peptide_data[np.invert(merged_flag_list)]
        filtered_df.to_csv(unknown_pep_file, sep=',', mode='w', header=True,
                           line_terminator='\n')


def main():
    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    csv_data = extract_csv_data("data/propep_g.csv")
    database_file = "data/sample_sprot.fasta"
    # if database_file.endswith(('fasta', 'fa')):
    #     with open(database_file, "r") as database:
    #         flags = search_peptide_db(csv_data, database)
    # else:
    #     with gzip.open(database_file, "rt") as database:
    #         flags = search_peptide_db(csv_data, database)
    pool = mp.Pool(processes=4)
    results = pool.map(search_peptide_db, [(csv_data, database_file, n) for n in range(1, 4 + 1)])
    pool.close()
    pool.join()
    # print(len(results))
    write_unknown_peptide_data(csv_data, results)
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main()
