#!/usr/bin/python3
"""A simple pipeline to filter peptide lists and preps data  for future use"""
import argparse
import datetime
import multiprocessing as mp
import os
import sys

import csv_dataframe
import peaks_peptide_comparison
import peptide_frequency
import peptide_venn
import unknown_peptide_seeker


def make_directories(dir_name):
    try:
        os.makedirs("output/")
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}".format(dir_name))
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}/comparison_output".format(dir_name))
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}/comparison_graphs".format(dir_name))
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}/peptide_count".format(dir_name))
    except FileExistsError:
        pass
    try:
        os.makedirs("output/{}/unknown_peptides".format(dir_name))
    except FileExistsError:
        pass


def compare_samples(left, right, name, left_name, right_name, prefix):
    try:
        print("*** Comparing {} and {} ***".format(left_name, right_name))
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

        peaks_peptide_comparison.find_distinct_peptides(left, right, prefix, left_name, right_name, name)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


def find_unknowns(dataframe, database, name, directory):
    try:
        print("*** finding {} peptides not appearing in database***".format(name))
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

        cpu = os.cpu_count()
        pool = mp.Pool(processes=cpu)
        results = pool.map(unknown_peptide_seeker.search_peptide_db,
                           [(dataframe, database, cpu, i) for i in range(1, cpu + 1)])
        pool.close()
        pool.join()

        merged_flags = unknown_peptide_seeker.merge_flags(results)
        unknown_peptide_seeker.write_unknown_peptide_data(dataframe, merged_flags, name, directory)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


def compare_distinct_unknown(directory, data_name):
    left = "output/{}/comparison_output/all_distinct_{}.csv".format(directory, data_name)
    left_data = csv_dataframe.extract_csv_data(left, drop_dupes=True)
    right = "output/{}/unknown_peptides/{}_unknown_peptides.csv".format(directory, data_name)
    right_data = csv_dataframe.extract_csv_data(right, drop_dupes=True)
    prefix = data_name
    left_name = "distinct"
    right_name = "unknown"
    compare_samples(left_data, right_data, directory, left_name, right_name, prefix)


def create_graph(left, right, left_name, right_name, directory):
    try:
        print("*** Creating Venn diagram ***")
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

        peptide_venn.create_venn_diagrams(left, right, left_name, right_name, directory)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files")
        print(e)
        sys.exit(2)


def create_sub_graph(directory, data_name):
    left = "output/{}/comparison_output/all_distinct_{}.csv".format(directory, data_name)
    left_data = csv_dataframe.extract_csv_data(left, drop_dupes=True)
    right = "output/{}/comparison_output/{}_distinct_distinct.csv".format(directory, data_name)
    right_data = csv_dataframe.extract_csv_data(right, drop_dupes=True)
    left_name = "distinct {}".format(data_name)
    right_name = "distinct database {}".format(data_name)
    create_graph(left_data, right_data, left_name, right_name, directory)


def count_peptides(left, right, left_name, right_name, directory):
    try:
        print("*** Counting peptides***")
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        peptides = peptide_frequency.create_peptide_list(left, right)
        peptide_frequency.create_counter_dataframe(left, left_name, directory, peptides)
        peptide_frequency.create_counter_dataframe(right, right_name, directory, peptides)

        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


def main(argv):
    print(' '.join(argv))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', '--left', action='store', dest="left",
                        help="Specify the .txt file containing the first group of peptide .csv file paths"
                             "containing paths to them", required=True)
    parser.add_argument('-r', '--right', action='store', dest="right",
                        help="Specify the .txt file containing the second group of peptide .csv file paths"
                             "containing paths to them", required=True)
    parser.add_argument('-d', '--database', action='store', dest="database", required=True,
                        help="Specify the directory of the protein database file")
    parser.add_argument('-n', '--name', action='store', dest="name", default="sample",
                        help="Provide a name for the pipeline run")
    parser.add_argument('--left_name', action='store', dest="left_name", default="left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action='store', dest="right_name", default="right",
                        help="Name the right sample")
    args = parser.parse_args()

    try:
        print("Pipeline started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        make_directories(args.name)
        left_data = csv_dataframe.join_dataframes(args.left)
        right_data = csv_dataframe.join_dataframes(args.right)

        compare_samples(left_data, right_data, args.name, args.left_name, args.right_name, "all")

        find_unknowns(left_data, args.database, args.left_name, args.name)
        find_unknowns(right_data, args.database, args.right_name, args.name)

        compare_distinct_unknown(args.name, args.left_name)
        compare_distinct_unknown(args.name, args.right_name)

        create_graph(left_data, right_data, args.left_name, args.right_name, args.name)
        create_sub_graph(args.name, args.left_name)
        create_sub_graph(args.name, args.right_name)

        count_peptides(args.left, args.right, args.left_name, args.right_name, args.name)

        print("Pipeline finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)