#!/usr/bin/python3
"""
doc
"""
import argparse
import datetime
import os
import sys

import pandas as pd
import scipy.stats as stats

import csv_dataframe


def join_dataframes(data):
    joined_dataframe = pd.DataFrame()
    with open(data, "r") as file_list:
        for file in file_list:
            csv_data = csv_dataframe.extract_csv_data(file.strip())
            joined_dataframe = joined_dataframe.append(csv_data[['Peptide']], ignore_index = True)
    return joined_dataframe


def count_peptide_frequency(peptide_data, column_name):
    count = pd.DataFrame(peptide_data['Peptide'].value_counts().reset_index())
    count.columns = ['Peptide', column_name]
    return count


def create_counter_dataframe(lefts, rights):
    all_pep = lefts.append(rights, ignore_index=True).drop_duplicates(subset=['Peptide'], keep='first') \
        .reset_index(drop=True)
    left_count = count_peptide_frequency(lefts, 'left_count')
    right_count = count_peptide_frequency(rights, 'right_count')

    merged = pd.merge(all_pep, left_count, on='Peptide', how='outer').fillna(0, downcast='infer')
    merged = pd.merge(merged, right_count, on='Peptide', how='outer').fillna(0, downcast='infer')

    with open("output/test1_rank.csv", "w+") as output_file:
        merged.to_csv(output_file, sep=',', mode='w', line_terminator='\n')

    return merged


def mann_whitney_u_test(merged_data):
    u_statistic, p_value = stats.mannwhitneyu(merged_data.left_count, merged_data.right_count)
    print('U-Statistic: ', u_statistic)
    print('p-value: ', p_value)


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
        left_data = join_dataframes(args.left)
        right_data = join_dataframes(args.right)
        merged_data = create_counter_dataframe(left_data, right_data)
        mann_whitney_u_test(merged_data)


if __name__ == '__main__':
    main(sys.argv)
