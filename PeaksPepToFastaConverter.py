#!/usr/bin/python3

import sys
import datetime
import pandas
import getopt
import os


def find_pep_from_accession(accession):
    """retrieves pep sequence from the full pep fasta file"""
    output = "distinct_pep_fasta/{}_distinct_pep.fasta".format("test")
    with open(output, "w+") as new_fasta:
        print("Hello, World!")


def extract_comparison_csv_data():
    """Reads the comparison csv data as dataframe"""
    input_file = "comparison_output/test1_distinct_gm.csv"
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    unique_accessions = csv_data[['Protein Accession']].drop_duplicates(keep='first')

    return unique_accessions


accessions = extract_comparison_csv_data()
