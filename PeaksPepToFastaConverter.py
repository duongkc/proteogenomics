#!/usr/bin/python3

import sys
import datetime
import pandas
import getopt
import os


def find_pep_from_accession(old_fasta, accession):
    """retrieves pep sequence from the full pep fasta file"""
    sequence = ""
    for line in old_fasta:
        if line.startswith(">") and accession in line:
            sequence = ""
            continue
        elif line.startswith(">") and accession not in line:
            break
        sequence += line
    return sequence


def extract_comparison_csv_data():
    """Reads the comparison csv data as dataframe"""
    input_file = "comparison_output/fake_sample.csv"
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    unique_accessions = csv_data[['Protein Accession']].drop_duplicates(keep='first')

    return unique_accessions


def main():
    try:
        os.makedirs("distinct_pep_fasta")
    except FileExistsError:
        pass

    accessions = extract_comparison_csv_data()
    output = "distinct_pep_fasta/{}_distinct_pep.fasta".format("test")
    old_pep_fasta = "output/Trinity.fasta.genemark.pep"
    with open(output, "w+") as new_fasta, open(old_pep_fasta, "r") as old_fasta:
        for i, row in accessions.iterrows():
            accession = row['Protein Accession']
            accession_seq = find_pep_from_accession(old_fasta, accession)
            new_fasta.write(">" + accession + "\n")
            new_fasta.write(accession_seq)


if __name__ == '__main__':
    main()


