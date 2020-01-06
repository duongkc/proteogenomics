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

import numpy as np
import pandas
from matplotlib import pyplot as plt


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
    csv_data = csv_data.drop_duplicates(subset=['Protein Accession', 'Peptide'], keep='first')
    return csv_data


def create_venn_diagrams(decoy_transdecoder, decoy_genemark, transdecoder, genemark, prefix):
    """Creates venn diagrams from the unique peptide lists of each file"""

    full_length_decoy_transdecoder = len(decoy_transdecoder.index)
    full_length_decoy_genemark = len(decoy_genemark.index)
    full_length_transdecoder = len(transdecoder.index)
    full_length_genemark = len(genemark.index)

    unique_length_decoy_transdecoder = len(
        decoy_transdecoder.drop_duplicates(subset=['Protein Accession'], keep='first'))
    unique_length_decoy_genemark = len(decoy_genemark.drop_duplicates(subset=['Protein Accession'], keep='first'))
    unique_length_transdecoder = len(transdecoder.drop_duplicates(subset=['Protein Accession'], keep='first'))
    unique_length_genemark = len(genemark.drop_duplicates(subset=['Protein Accession'], keep='first'))

    tdd_merged = pandas.merge(decoy_transdecoder, decoy_genemark, on='Peptide', how='left', indicator=True) \
        .query("_merge == 'left_only'")
    gmd_merged = pandas.merge(decoy_transdecoder, decoy_genemark, on='Peptide', how='right', indicator=True) \
        .query("_merge == 'right_only'")
    td_merged = pandas.merge(transdecoder, genemark, on='Peptide', how='left', indicator=True) \
        .query("_merge == 'left_only'")
    gm_merged = pandas.merge(transdecoder, genemark, on='Peptide', how='right', indicator=True) \
        .query("_merge == 'right_only'")

    unique_merged_decoy_transdecoder = len(tdd_merged.drop_duplicates(subset=['Protein Accession_x'], keep='first'))
    unique_merged_decoy_genemark = len(gmd_merged.drop_duplicates(subset=['Protein Accession_y'], keep='first'))
    unique_merged_transdecoder = len(td_merged.drop_duplicates(subset=['Protein Accession_x'], keep='first'))
    unique_merged_genemark = len(gm_merged.drop_duplicates(subset=['Protein Accession_y'], keep='first'))

    full_lengths = [full_length_transdecoder, full_length_genemark,
                    full_length_decoy_transdecoder, full_length_decoy_genemark]
    unique_lengths = [unique_length_transdecoder, unique_length_genemark,
                      unique_length_decoy_transdecoder, unique_length_decoy_genemark]
    merged_gms = [len(gm_merged), len(gmd_merged)]
    merged_tds = [len(td_merged), len(tdd_merged)]

    unique_merged_tds = [unique_merged_transdecoder, unique_merged_decoy_transdecoder]
    unique_merged_gms = [unique_merged_genemark, unique_merged_decoy_genemark]
    names = ['Transdecoder', 'GenemarkS-T', 'Transdecoder+decoy', 'GenemarkS-T+decoy']
    fig, axes = plt.subplots(2, 2, figsize=(15, 8))
    axes[0][0].bar(names, full_lengths, color=('#67C3FF', '#FF8267'))
    axes[0][1].bar(names, unique_lengths, color=('#67C3FF', '#FF8267'))

    index = np.arange(2)
    rects1 = axes[1][0].bar(index, merged_tds, 0.35, color='#67C3FF', label='Transdecoder')
    rects2 = axes[1][0].bar(index + 0.35, merged_gms, 0.35, color='#FF8267', label='GenemarkS-T')

    axes[1][0].set_xticks(index + 0.175)
    axes[1][0].set_xticklabels(['Without decoy', 'With decoy'])
    axes[1][0].legend()

    rects3 = axes[1][1].bar(index, unique_merged_tds, 0.35, color='#67C3FF', label='Transdecoder')
    rects4 = axes[1][1].bar(index + 0.35, unique_merged_gms, 0.35, color='#FF8267', label='GenemarkS-T')

    axes[1][1].set_xticks(index + 0.175)
    axes[1][1].set_xticklabels(['Without decoy', 'With decoy'])
    axes[1][1].legend()

    axes[0][0].title.set_text('Protein-peptide matches')
    axes[0][1].title.set_text('Unique proteins in the matches')
    axes[1][0].title.set_text('Protein-peptide matches unique to database')
    axes[1][1].title.set_text('Unique proteins in matches unique to database')

    plt.suptitle('Sample {} Protein-Peptide matches'.format(prefix))
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.tight_layout()
    plt.savefig('output/comparison_graphs/sample_{}_bar.png'.format(prefix))


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
        os.makedirs("output/comparison_graphs")
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
