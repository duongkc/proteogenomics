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
    with open(data, "r") as file_list:
        for file in file_list:
            csv_data = csv_dataframe.extract_csv_data(file)



def main(argv):
    print(' '.join(argv))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--batch', action="store_true")
    parser.add_argument('-l', '--left', action="store", dest="left",
                        help="Specify the directory of the first protein-peptides.csv file", required=True)
    parser.add_argument('-r', '--right', action="store", dest="right",
                        help="Specify the directory of the second protein-peptides.csv file", required=True)
    parser.add_argument('-p', '--prefix', action="store", dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    parser.add_argument('--left_name', action="store", dest="left_name", default="left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action="store", dest="right_name", default="right",
                        help="Name the right sample")
    args = parser.parse_args()


if __name__ == '__main__':
    main(sys.argv)
