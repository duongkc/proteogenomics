#!/usr/bin/python3
"""
usage: venn_test.py -g <genemark csv file> -t <transdecoder csv file>
                  -h <decoy genemark file> -u <decoy transdecoder file> -p <output prefix>
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
    plt.savefig('comparison_graphs/sample_{}.png'.format(prefix))


def main(argv):
    print(' '.join(argv))
    real_genemark_file = ''
    real_trans_file = ''
    decoy_genemark_file = ''
    decoy_trans_file = ''
    output_prefix = ''

    try:
        opts, args = getopt.getopt(argv[1:], 'g:t:h:u:p:', ['genemark=', 'transdecoder=',
                                                            'genemark_decoy=', 'transdecoder_decoy=', 'prefix='])
    except (getopt.GetoptError, FileNotFoundError) as e:
        print("usage: venn_test.py -g <genemark csv file> -t <transdecoder csv file> "
              "-h <decoy genemark file> -u <decoy transdecoder file> -p <output prefix>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-g', '--genemark'):
            real_genemark_file = arg
        elif opt in ('-t', '--transdecoder'):
            real_trans_file = arg
        elif opt in ('-h', '--genemark_decoy'):
            decoy_genemark_file = arg
        elif opt in ('-u', '--transdecoder_decoy'):
            decoy_trans_file = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print("usage: venn_test.py -g <genemark csv file> -t <transdecoder csv file> "
                  "-h <decoy genemark file> -u <decoy transdecoder file> -p <output prefix>")
            sys.exit(2)

    try:
        os.makedirs("comparison_graphs")
    except FileExistsError:
        pass

    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    decoy_trans_data = extract_csv_data(decoy_trans_file)
    decoy_genemark_data = extract_csv_data(decoy_genemark_file)
    real_trans_data = extract_csv_data(real_trans_file)
    real_genemark_data = extract_csv_data(real_genemark_file)

    create_venn_diagrams(decoy_trans_data, decoy_genemark_data, real_trans_data, real_genemark_data, output_prefix)
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main(sys.argv)
