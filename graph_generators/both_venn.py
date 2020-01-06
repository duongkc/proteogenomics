#!/usr/bin/python3
"""
Creates venn diagrams in one image, but only for comparing GenemarkS-T with Transdecoder csv files, and not
the decoy files.

Usage:

"""

import datetime
import getopt
import os
import re
import sys

import pandas
from matplotlib import pyplot as plt
from matplotlib_venn import venn2


def clean_peptide_col(peptide_column):
    """Cleans up the peptide column by removing unnecessary information and returns the peptide"""
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file):
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    csv_data = csv_data[['Protein Accession', 'Peptide']].drop_duplicates(subset=['Peptide'], keep='first')
    return csv_data


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
        trans_data = extract_csv_data(trans_file)
        genemark_data = extract_csv_data(genemark_file)

        create_venn_diagrams(genemark_data, trans_data, output_prefix)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(e)
        print(__doc__)

        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)