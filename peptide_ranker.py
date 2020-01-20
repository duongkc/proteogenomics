#!/usr/bin/python3
"""
doc
"""
import argparse
import datetime
import os
import sys

import pandas

import csv_dataframe


def join_dataframes(data):
    joined_dataframe = pandas.DataFrame()
    with open(data, "r") as file_list:
        for file in file_list:
            csv_data = csv_dataframe.extract_csv_data(file.strip())
            joined_dataframe = joined_dataframe.append(csv_data[['Peptide']], ignore_index = True)
    return joined_dataframe


def main(argv):
    print(' '.join(argv))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--batch', action='store_true', dest="batch",
                        help="Include this argument if the left and right csv files are specified in a batch txt file")
    parser.add_argument('-l', '--left', action='store', dest="left",
                        help="Specify the directory of the first protein-peptides.csv file or batch .txt"
                             "containing paths to them", required=True)
    parser.add_argument('-r', '--right', action='store', dest="right",
                        help="Specify the directory of the second protein-peptides.csv file or batch .txt"
                             "containing paths to them", required=True)
    parser.add_argument('-p', '--prefix', action='store', dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    parser.add_argument('--left_name', action='store', dest="left_name", default="left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action='store', dest="right_name", default="right",
                        help="Name the right sample")
    args = parser.parse_args()

    if args.batch:
        lefts = join_dataframes(args.left)
        rights = join_dataframes(args.right)
        print(lefts['Peptide'].value_counts())
        # print(type(rights['Peptide'].value_counts()))
        with open("output/test_rank.csv", "w+") as output_file:
            lefts.to_csv(output_file, sep=',', mode='w', line_terminator='\n')


if __name__ == '__main__':
    main(sys.argv)
