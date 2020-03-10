#!/usr/bin/python3
"""
A module that checks how many peptides from the PEAKS protein-peptide match results can be found in
existing protein sequence databases and filters the unknown peptides to a separate file.
"""
import argparse
import datetime
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


def write_unknown_peptide_data(peptide_data, merged_flag_list, prefix, directory):
    """Writes filtered dataframe to a new csv file"""
    output = "output/{}/unknown_peptides/{}_unknown_peptides.csv".format(directory, prefix)

    with open(output, "w") as unknown_pep_file:
        filtered_df = peptide_data[merged_flag_list]
        filtered_df.to_csv(unknown_pep_file, sep=',', mode='w', header=True,
                           line_terminator='\n')


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', '--csv', action='store', dest="csv", required=True,
                        help="Specify the .txt file containing the peptide .csv file paths")
    parser.add_argument('-d', '--database', action='store', dest="database", required=True,
                        help="Specify the directory of the protein database file")
    parser.add_argument('-o', '--outdir', action='store', dest='outdir', default="peptides",
                        help="Provide an output directory name, i.e. 'output/<NAME>/unknown_peptides/'")
    parser.add_argument('-p', '--prefix', action='store', dest="prefix", default="sample",
                        help="Provide a prefix for the output file: <PREFIX>_unknown_peptides.csv")

    args = parser.parse_args()

    try:
        os.makedirs("output/{}/unknown_peptides".format(args.outdir))
    except FileExistsError:
        pass

    try:
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

        csv_data = csv_dataframe.join_dataframes(args.csv)

        cpu = os.cpu_count()
        pool = mp.Pool(processes=cpu)
        results = pool.map(search_peptide_db, [(csv_data, args.database, cpu, i) for i in range(1, cpu + 1)])
        pool.close()
        pool.join()

        merged_flags = merge_flags(results)
        write_unknown_peptide_data(csv_data, merged_flags, args.prefix, args.outdir)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
