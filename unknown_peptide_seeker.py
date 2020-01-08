#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.

Usage: unknown_peptide_seeker.py -c <PEAKS csv file> -d <protein database fasta> -p <output prefix>
"""
import datetime
import getopt
import multiprocessing as mp
import os
import sys

import numpy as np
from Bio import SeqIO

import csv_dataframe


def search_peptide_db(arguments):
    """Checks for presence of peptides in protein database and notes them down in a boolean list"""
    peptide_data, database_file, n, offset = arguments

    with open(database_file, "r") as database:
        count = 1
        flag_list = [1] * len(peptide_data.index)
        for record in SeqIO.parse(database, "fasta"):
            if (count - offset) % n == 0:
                for i, row in peptide_data.iterrows():
                    peptide = row['Peptide']
                    if peptide in record.seq:
                        flag_list[i] = 0
            count += 1
    return flag_list


def merge_flags(flags):
    """Merges the boolean lists of the separate processes into 1 list"""
    for i in range(len(flags)):
        flags[i] = np.array(flags[i], dtype=bool)

    merged_flag_list = flags[0]
    for i in range(1, len(flags)):
        merged_flag_list = merged_flag_list & flags[i]
    return merged_flag_list


def write_unknown_peptide_data(peptide_data, merged_flag_list, prefix):
    """Writes filtered dataframe to a new csv file"""
    output = "output/unknown_peptides/{}_unknown_peptides.csv".format(prefix)

    with open(output, "w") as unknown_pep_file:
        filtered_df = peptide_data[merged_flag_list]
        filtered_df.to_csv(unknown_pep_file, sep=',', mode='w', header=True,
                           line_terminator='\n')


def main(argv):
    print(' '.join(argv))
    print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

    csv_file = ''
    database_file = ''
    output_prefix = ''
    cpu = os.cpu_count()

    try:
        opts, args = getopt.getopt(argv[1:], 'c:d:p:', ['csv=', 'database=', 'prefix='])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-c', '--csv'):
            csv_file = arg
        elif opt in ('-d', '--database'):
            database_file = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print(__doc__)
            sys.exit(2)

    try:
        os.makedirs("output/unknown_peptides")
    except FileExistsError:
        pass

    try:
        csv_data = csv_dataframe.extract_csv_data(csv_file)
        # if database_file.endswith(('fasta', 'fa')):
        #     with open(database_file, "r") as database:
        #         flags = search_peptide_db(csv_data, database)
        # else:
        #     with gzip.open(database_file, "rt") as database:
        #         flags = search_peptide_db(csv_data, database)
        pool = mp.Pool(processes=cpu)
        results = pool.map(search_peptide_db, [(csv_data, database_file, cpu, i) for i in range(1, cpu + 1)])
        pool.close()
        pool.join()

        merged_flags = merge_flags(results)
        write_unknown_peptide_data(csv_data, merged_flags, output_prefix)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError:
        print(__doc__)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
