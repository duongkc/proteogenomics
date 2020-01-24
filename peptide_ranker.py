#!/usr/bin/python3
"""
A python module that counts and compares peptide frequency between two groups. Results can be found in a separate
output file.
A Mann-Whitney U test will also be performed with these results, with the results shown in the terminal.
If specified, a Wilcoxon signed-rank test will be performed instead.
"""
import argparse
import datetime
import os
import sys

import pandas as pd
import scipy.stats as stats
import numpy as np

import csv_dataframe


def join_dataframes(data):
    joined_dataframe = pd.DataFrame()
    with open(data, "r") as file_list:
        for file in file_list:
            csv_data = csv_dataframe.extract_csv_data(file.strip(), drop_dupes=True)
            joined_dataframe = joined_dataframe.append(csv_data[['Peptide']], ignore_index=True)
    return joined_dataframe


def create_peptide_list(left_file, right_file):
    """Creates a list of all peptides as a DataFrame column"""
    joined_left = join_dataframes(left_file)
    joined_right = join_dataframes(right_file)
    all_peptides = joined_left.append(joined_right, ignore_index=True) \
        .drop_duplicates(subset=['Peptide'], keep='first').reset_index(drop=True)
    with open("output/all_peptides_gm.csv", "w+") as output:
        all_peptides.to_csv(output, sep=',', mode='w', line_terminator='\n')


def count_peptide_frequency(peptide_data, column_name):
    count = pd.DataFrame(peptide_data['Peptide'].value_counts().reset_index())
    count.columns = ['Peptide', column_name]
    return count


def parts_per_million(data):
    data_sum = data.sum()
    ppm = data / data_sum * 1000000
    return ppm


def create_counter_dataframe(files, group_name, prefix):
    output_file = "output/{}_peptide_frequency_{}.csv".format(prefix, group_name)
    all_peptides = pd.read_csv("output/all_peptides_gm.csv", header='infer', delimiter=',', index_col=0)
    with open(files, "r") as file_list:
        for num, file in enumerate(file_list):
            file_data = csv_dataframe.extract_csv_data(file.strip(), drop_dupes=False)
            counter_column = count_peptide_frequency(file_data, "{}{}".format(group_name[0].capitalize(), num + 1))
            all_peptides = pd.merge(all_peptides, counter_column, on='Peptide', how='outer')

    all_peptides = all_peptides.fillna(0, downcast='infer')
    for column in all_peptides.columns[1:]:
        all_peptides[column] = parts_per_million(all_peptides[column])
    # merged['abs'] = np.abs(merged['left_count'] - merged['right_count'])
    # merged = merged.sort_values(by=['abs'], ascending=False)
    #
    with open(output_file, "w+") as output_file:
        all_peptides.to_csv(output_file, sep=',', mode='w', line_terminator='\n')
    return all_peptides


def mann_whitney_u_test(left_data, right_data):
    peptides = left_data[['Peptide']].copy()
    for i, row in left_data.iterrows():
        left = left_data.iloc[i, 1:]
        right = right_data.iloc[i, 1:]
        if all(sample == left[0] for sample in left):
            if all(sample == right[0] for sample in right):
                continue
        u_statistic, p_value = stats.mannwhitneyu(left, right, alternative='two-sided')
        peptides.at[i, 'p-value'] = p_value
        peptides.at[i, 'u-statistic'] = u_statistic
    peptides = peptides.sort_values(by=['p-value'], ascending=True)
    with open("output/mann_peptides.csv", "w+") as output:
        peptides.to_csv(output, sep=',', mode='w', line_terminator='\n')


def wilcoxon_test(merged_data):
    w_statistic, p_value = stats.wilcoxon(merged_data.left_count, merged_data.right_count, alternative='two-sided')
    print('W-Statistic: ', w_statistic)
    print('p-value: ', p_value)


def main(argv):
    print(' '.join(argv))
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-l', '--left', action='store', dest="left",
                        help="Specify the .txt file containing the first group of peptide .csv file paths"
                             "containing paths to them", required=True)
    parser.add_argument('-r', '--right', action='store', dest="right",
                        help="Specify the .txt file containing the second group of peptide .csv file paths"
                             "containing paths to them", required=True)
    parser.add_argument('-p', '--prefix', action='store', dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    parser.add_argument('--left_name', action='store', dest="left_name", default="left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action='store', dest="right_name", default="right",
                        help="Name the right sample")
    parser.add_argument('-w', '--wilcoxon', action='store_true', dest="wilcoxon",
                        help="Will perform a Wilcoxon signed-rank test instead of a Mann-Whitney U test if this"
                             "argument is provided.")
    args = parser.parse_args()

    try:
        print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        left_data = create_counter_dataframe(args.left, args.left_name, args.prefix)
        right_data = create_counter_dataframe(args.right, args.right_name, args.prefix)
        mann_whitney_u_test(left_data, right_data)
        # if args.batch:
        #     left_data = join_dataframes(args.left)
        #     right_data = join_dataframes(args.right)
        #     merged_data = create_counter_dataframe(left_data, right_data, args.prefix)
        #     if args.wilcoxon:
        #         wilcoxon_test(merged_data)
        #     else:
        #         mann_whitney_u_test(merged_data)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
