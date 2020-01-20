#!/usr/bin/python3
"""
Creates venn diagrams in one image, by comparing PEAKS protein-peptide.csv files or any csv file using that format.
"""

import argparse
import datetime
import os
import sys

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

import csv_dataframe


def create_venn_diagrams(left, right, suffix, left_name, right_name):
    """Creates venn diagrams from the unique peptide lists of each file"""

    set_left = set(left.Peptide)
    set_right = set(right.Peptide)

    fig, axes = plt.subplots(nrows=1, ncols=2)
    v1 = venn2([set_left, set_right], set_labels=(left_name, right_name), ax=axes[0])

    total_v2 = len(set_left.union(set_right))
    v2 = venn2([set_left, set_right], set_labels=(left_name, right_name), ax=axes[1],
               subset_label_formatter=lambda x: f"{(x/total_v2):.2%}")

    plt.suptitle('Comparing {} peptide matches'.format(suffix))
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.tight_layout()
    plt.savefig('output/comparison_graphs/venn_{}.png'.format(suffix))


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-l', '--left', action='store', dest="left",
                        help="Specify directory of the first sample .csv", required=True)
    parser.add_argument('-r', '--right', action='store', dest="right",
                        help="Specify directory of the second sample .csv", required=True)
    parser.add_argument('-s', '--suffix', action='store', dest="suffix", default="peptides",
                        help="Give the output file a custom name: venn_<suffix>.png")
    parser.add_argument('--left_name', action='store', dest="left_name", default="sample left",
                        help="Name the left sample")
    parser.add_argument('--right_name', action='store', dest="right_name", default="sample right",
                        help="Name the right sample")
    args = parser.parse_args()

    try:
        os.makedirs("output/comparison_graphs")
    except FileExistsError:
        pass

    try:
        print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        left_data = csv_dataframe.extract_csv_data(args.left)
        right_data = csv_dataframe.extract_csv_data(args.right)

        create_venn_diagrams(left_data, right_data, args.suffix, args.left_name, args.right_name)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files")

        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
