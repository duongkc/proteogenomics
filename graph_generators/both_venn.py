#!/usr/bin/python3
"""
Creates venn diagrams in one image, but only for comparing GenemarkS-T with Transdecoder csv files, and not
the decoy files.

Usage:

"""

import datetime
import getopt
import os
import sys

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

import csv_dataframe


def create_venn_diagrams(genemark, transdecoder, prefix):
    """Creates venn diagrams from the unique peptide lists of each file"""

    set_td = set(transdecoder.Peptide)
    set_gm = set(genemark.Peptide)

    fig, axes = plt.subplots(nrows=1, ncols=2)
    v1 = venn2([set_td, set_gm], set_labels=('Sample 1', 'Sample 2'), ax=axes[0])

    total_v2 = len(set_td.union(set_gm))
    v2 = venn2([set_td, set_gm], set_labels=('Sample 1', 'Sample 2'), ax=axes[1],
               subset_label_formatter=lambda x: f"{(x/total_v2):.2%}")

    plt.suptitle('Comparing Sample {} peptide matches'.format(prefix))
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.tight_layout()
    plt.savefig('output/comparison_graphs/venn_sample_{}.png'.format(prefix))


def main(argv):
    print(' '.join(argv))
    genemark_file = ''
    trans_file = ''
    output_prefix = ''

    try:
        opts, args = getopt.getopt(argv[1:], 'g:t:p:', ['genemark=', 'transdecoder=', 'prefix='])
    except (getopt.GetoptError, FileNotFoundError) as e:
        print(__doc__)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-g', '--genemark'):
            genemark_file = arg
        elif opt in ('-t', '--transdecoder'):
            trans_file = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print(__doc__)
            sys.exit(2)

    try:
        os.makedirs("output/comparison_graphs")
    except FileExistsError:
        pass

    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    try:
        trans_data = csv_dataframe.extract_csv_data(trans_file)
        genemark_data = csv_dataframe.extract_csv_data(genemark_file)

        create_venn_diagrams(genemark_data, trans_data, output_prefix)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(e)
        print(__doc__)

        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)