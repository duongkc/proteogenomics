#!/usr/bin/python3

import sys
import datetime
import pandas
import getopt
import os

"""
usage: PeaksPepToFastaConverter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>
"""


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


def extract_comparison_csv_data(input_file):
    """Reads the comparison csv data as dataframe"""
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    unique_accessions = csv_data[['Protein Accession']].drop_duplicates(keep='first')

    return unique_accessions


def main(argv):
    print(' '.join(argv))
    print("started at: " + datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))
    input_csv = ''
    old_pep_fasta = ''
    output_prefix = ''

    try:
        opts, args = getopt.getopt(argv[1:], 'c:f:p:', ['csv=', 'fasta=', 'prefix='])
    except getopt.GetoptError:
        print("usage: PeaksPepToFastaConverter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-c', '--csv'):
            input_csv = arg
        elif opt in ('-f', '--fasta'):
            old_pep_fasta = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print("usage: PeaksPepToFastaConverter.py -c <distinct peptide csv> -f <peptide fasta> -p <output prefix>")
            sys.exit(2)

    try:
        os.makedirs("distinct_pep_fasta")
    except FileExistsError:
        pass

    accessions = extract_comparison_csv_data(input_csv)
    output = "distinct_pep_fasta/{}_distinct_pep.fasta".format(output_prefix)
    with open(output, "w+") as new_fasta, open(old_pep_fasta, "r") as old_fasta:
        for i, row in accessions.iterrows():
            accession = row['Protein Accession']
            accession_seq = find_pep_from_accession(old_fasta, accession)
            new_fasta.write(">" + accession + "\n")
            new_fasta.write(accession_seq)
    print("finished at: " + datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main(sys.argv)


