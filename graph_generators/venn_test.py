#!/usr/bin/python3
"""
usage: venn_test.py -g <genemark csv file> -t <transdecoder csv file>
                  -h <decoy genemark file> -u <decoy transdecoder file> -p <output prefix>
"""

import argparse
import datetime
import os
import sys

from matplotlib import pyplot as plt
from matplotlib_venn import venn2

import csv_dataframe


def create_venn_diagrams(decoy_transdecoder, decoy_genemark, transdecoder, genemark, prefix):
    """Creates venn diagrams from the unique peptide lists of each file"""

    set_td_decoy = set(decoy_genemark.Peptide)
    set_gm_decoy = set(decoy_transdecoder.Peptide)
    set_td = set(transdecoder.Peptide)
    set_gm = set(genemark.Peptide)

    fig, axes = plt.subplots(nrows=2, ncols=2)
    # total_v1 = len(list_td.union(list_gm))
    # v1 = venn2([list_td, list_gm], set_labels=('Transdecoder', 'GenemarkS-T'), ax=axes[0][0],
    #            subset_label_formatter=lambda x: f"{(x/total_v1):1.0%}")
    v1 = venn2([set_td, set_gm], set_labels=('Transdecoder', 'GenemarkS-T'), ax=axes[0][0])
    v2 = venn2([set_td_decoy, set_gm_decoy], set_labels=('Transdecoder decoy', 'GenemarkS-T decoy'), ax=axes[0][1])
    v3 = venn2([set_td, set_td_decoy], set_labels=('Transdecoder', 'Transdecoder decoy'), ax=axes[1][0])
    v4 = venn2([set_gm, set_gm_decoy], set_labels=('GenemarkS-T', 'GenemarkS-T decoy'), ax=axes[1][1])

    plt.suptitle('Sample {} peptide matches'.format(prefix))
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.tight_layout()
    plt.savefig('output/comparison_graphs/sample_{}.png'.format(prefix))


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-g', '--genemark', action="store", dest="genemark", required=True,
                        help="Specify the directory of the GenemarkS-T protein-peptides.csv file")
    parser.add_argument('-t', '--transdecoder', action="store", dest="transdecoder", required=True,
                        help="Specify the directory of the Transdecoder protein-peptides.csv file")
    parser.add_argument('-h', '--genemark_decoy', action="store", dest="genemark_decoy", required=True,
                        help="Specify the directory of the decoy GenemarkS-T protein-peptides.csv file")
    parser.add_argument('-u', '--transdecoder_decoy', action="store", dest="transdecoder_decoy", required=True,
                        help="Specify the directory of the decoy Transdecoder protein-peptides.csv file")
    parser.add_argument('-p', '--prefix', action="store", dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    args = parser.parse_args()

    try:
        os.makedirs("output/comparison_graphs")
    except FileExistsError:
        pass

    try:
        print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        decoy_trans_data = csv_dataframe.extract_csv_data(args.transdecoder_decoy)
        decoy_genemark_data = csv_dataframe.extract_csv_data(args.genemark_decoy)
        real_trans_data = csv_dataframe.extract_csv_data(args.transdecoder)
        real_genemark_data = csv_dataframe.extract_csv_data(args.genemark)

        create_venn_diagrams(decoy_trans_data, decoy_genemark_data, real_trans_data, real_genemark_data, args.prefix)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError:
        print(__doc__)
        print("Please provide valid files")
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
