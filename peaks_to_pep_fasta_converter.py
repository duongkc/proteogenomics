#!/usr/bin/python3
"""Small module for converting the distinct matched gene accessions back to peptide fastas
usage: peaks_to_pep_fasta_converter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>
"""


import datetime
import getopt
import os
import sys

import pandas
from Bio import SeqIO


def write_found_peptides_to_fasta(old_fasta, new_fasta, accessions):
    """Writes every peptide record that matches the accession list to a new fasta file"""
    for record in SeqIO.parse(old_fasta, 'fasta'):
        if record.id.split()[0] in accessions:
            SeqIO.write(record, new_fasta, "fasta")


def extract_comparison_csv_data(input_file):
    """Reads the comparison csv data as dataframe"""
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    unique_accessions = csv_data['Protein Accession'].drop_duplicates(keep='first').values.tolist()

    return unique_accessions


def main(argv):
    print(' '.join(argv))
    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    input_csv = ''
    old_pep_fasta = ''
    output_prefix = ''

    try:
        opts, args = getopt.getopt(argv[1:], 'c:f:p:', ['csv=', 'fasta=', 'prefix='])
    except getopt.GetoptError:
        print('usage: peaks_to_pep_fasta_converter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-c', '--csv'):
            input_csv = arg
        elif opt in ('-f', '--fasta'):
            old_pep_fasta = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print(
                'usage: peaks_to_pep_fasta_converter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>'
            )
            sys.exit(2)

    try:
        os.makedirs("distinct_pep_fasta")
    except FileExistsError:
        pass

    accessions = extract_comparison_csv_data(input_csv)
    output = "distinct_pep_fasta/{}_distinct_pep.fasta".format(output_prefix)
    with open(output, "w+") as new_fasta, open(old_pep_fasta, "r") as old_fasta:
        write_found_peptides_to_fasta(old_fasta, new_fasta, accessions)
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main(sys.argv)


