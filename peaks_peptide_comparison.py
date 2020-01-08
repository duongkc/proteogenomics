#!/usr/bin/python3
"""Short module to compare two PEAKS psm search output files to find overlap and distinction between the two in terms
   of matched peptides.

    usage: peaks_peptide_comparison.py -g <genemark csv file> -t <transdecoder csv file> -p <output prefix>
"""


import datetime
import getopt
import os
import sys

import pandas

import csv_dataframe


def find_distinct_peptides(transdecoder_data, genemark_data, prefix):
    """Filters the CSV files so only distinct peptides remain"""
    distinct_td_csv = "output/comparison_output/{}_distinct_td.csv".format(prefix)
    distinct_gm_csv = "output/comparison_output/{}_distinct_gm.csv".format(prefix)
    with open(distinct_td_csv, "w+") as distinct_transdecoder, \
            open(distinct_gm_csv, "w+") as distinct_genemark:
        td_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='left', indicator=True) \
            .query("_merge == 'left_only'")
        gm_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='right', indicator=True) \
            .query("_merge == 'right_only'")
        td_merged[['Peptide']] \
            .to_csv(distinct_transdecoder, sep=',', mode='w', index=False, header=['Peptide'],
                    line_terminator='\n')
        gm_merged[['Peptide']] \
            .to_csv(distinct_genemark, sep=',', mode='w', index=False, header=['Peptide'],
                    line_terminator='\n')


def main(argv):
    print(' '.join(argv))
    genemark_csv = ''
    transdecoder_csv = ''
    output_prefix = ''

    try:
        opts, args = getopt.getopt(argv[1:], 'g:t:p:', ['genemark=', 'transdecoder=', 'prefix='])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-g', '--genemark'):
            genemark_csv = arg
        elif opt in ('-t', '--transdecoder'):
            transdecoder_csv = arg
        elif opt in ('-p', '--prefix'):
            output_prefix = arg
        else:
            print(__doc__)
            sys.exit(2)

    try:
        os.makedirs("output/comparison_output")
    except FileExistsError:
        pass

    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    genemark_data = csv_dataframe.extract_csv_data(genemark_csv)
    transdecoder_data = csv_dataframe.extract_csv_data(transdecoder_csv)

    find_distinct_peptides(transdecoder_data, genemark_data, output_prefix)
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main(sys.argv)
