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
    for text in v1.set_labels:
        text.set_color('w')
    for text in v1.subset_labels:
        text.set_color('w')
        text.set_fontsize(16)
    v1.get_patch_by_id('10').set_color('#FFD54D')
    v1.get_patch_by_id('10').set_edgecolor('none')
    v1.get_patch_by_id('10').set_alpha(0.8)
    v1.get_patch_by_id('01').set_color('#DF3A42')
    v1.get_patch_by_id('01').set_edgecolor('none')
    v1.get_patch_by_id('01').set_alpha(0.8)
    v1.get_patch_by_id('11').set_color('#EF7D24')
    v1.get_patch_by_id('11').set_edgecolor('none')
    v1.get_patch_by_id('11').set_alpha(1)
    v1.get_label_by_id('01').set_y(0.05)
    v1.get_label_by_id('11').set_y(-0.05)

    total_v2 = len(set_left.union(set_right))
    v2 = venn2([set_left, set_right], set_labels=(left_name, right_name), ax=axes[1],
               subset_label_formatter=lambda x: f"{(x/total_v2):.1%}")
    for text in v2.set_labels:
        text.set_color('w')
    for text in v2.subset_labels:
        text.set_color('w')
        text.set_fontsize(16)
    v2.get_patch_by_id('10').set_color('#FFD54D')
    v2.get_patch_by_id('10').set_edgecolor('none')
    v2.get_patch_by_id('10').set_alpha(0.8)
    v2.get_patch_by_id('01').set_color('#DF3A42')
    v2.get_patch_by_id('01').set_edgecolor('none')
    v2.get_patch_by_id('01').set_alpha(0.8)
    v2.get_patch_by_id('11').set_color('#EF7D24')
    v2.get_patch_by_id('11').set_edgecolor('none')
    v2.get_patch_by_id('11').set_alpha(1)
    v2.get_label_by_id('01').set_y(0.05)
    v2.get_label_by_id('11').set_y(-0.05)

    plt.suptitle('Comparing {} peptide matches'.format(suffix), color='w')
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.tight_layout()
    plt.savefig('output/comparison_graphs/venn_{}.png'.format(suffix), transparent=True)


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
        left_data = csv_dataframe.extract_csv_data(args.left, drop_dupes=True)
        right_data = csv_dataframe.extract_csv_data(args.right, drop_dupes=True)

        create_venn_diagrams(left_data, right_data, args.suffix, args.left_name, args.right_name)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
