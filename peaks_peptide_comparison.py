#!/usr/bin/python3
"""Short module to compare two PEAKS psm protein-peptide output files to find overlap and distinction between the two
    in terms of matched peptides. Distinct peptides of each file will be stored in separate files.
"""
import argparse
import datetime
import os
import sys

import pandas

import csv_dataframe


def find_distinct_peptides(left_data, right_data, prefix, left_name, right_name, output_dir):
    """Filters the CSV files so only distinct peptides remain"""
    distinct_left_csv = "output/{}/comparison_output/{}_distinct_{}.csv".format(output_dir, prefix, left_name)
    distinct_right_csv = "output/{}/comparison_output/{}_distinct_{}.csv".format(output_dir, prefix, right_name)
    with open(distinct_left_csv, "w+") as distinct_left, \
            open(distinct_right_csv, "w+") as distinct_right:
        left_merged = pandas.merge(left_data, right_data, on='Peptide', how='left', indicator=True) \
            .query("_merge == 'left_only'")
        right_merged = pandas.merge(left_data, right_data, on='Peptide', how='right', indicator=True) \
            .query("_merge == 'right_only'")
        left_merged[['Peptide']] \
            .to_csv(distinct_left, sep=',', mode='w', index=False, header=['Peptide'],
                    line_terminator='\n')
        right_merged[['Peptide']] \
            .to_csv(distinct_right, sep=',', mode='w', index=False, header=['Peptide'],
                    line_terminator='\n')
    with open("output/{}/comparison_output/{}_common_peps.csv".format(output_dir, prefix), "w+") as commons:
        common_peps = pandas.merge(left_data, right_data, on='Peptide', how='outer', indicator=True) \
            .query("_merge == 'both'")
        common_peps[['Peptide']]\
            .to_csv(commons, sep=',', mode='w', index=False, header=['Peptide'], line_terminator='\n')


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', '--left', action='store', dest="left",
                        help="Specify the .txt file containing the first group of peptide .csv file paths",
                        required=True)
    parser.add_argument('-r', '--right', action='store', dest="right",
                        help="Specify the .txt file containing the second group of peptide .csv file paths",
                        required=True)
    parser.add_argument('--left_name', action='store', dest="left_name", default="left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action='store', dest="right_name", default="right",
                        help="Name the right sample")
    parser.add_argument('-o', '--outdir', action='store', dest='outdir', default="peptides",
                        help="Provide an output directory name, i.e. 'output/<NAME>/comparison_output/'")
    parser.add_argument('-p', '--prefix', action='store', dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    args = parser.parse_args()

    try:
        os.makedirs("output/{}".format(args.outdir))
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}/comparison_output".format(args.outdir))
    except FileExistsError:
        pass

    try:
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        left_data = csv_dataframe.join_dataframes(args.left, batch=False)
        right_data = csv_dataframe.join_dataframes(args.right, batch=False)

        find_distinct_peptides(left_data, right_data, args.prefix, args.left_name, args.right_name, args.outdir)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
