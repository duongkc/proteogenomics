#!/usr/bin/python3
"""Short module to compare two PEAKS psm search output files to find overlap and distinction between the two in terms
   of matched peptides.

    usage: peaks_peptide_comparison.py -g <genemark csv file> -t <transdecoder csv file> -p <output prefix>
"""
import argparse
import datetime
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

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-g', '--genemark', action="store", dest="genemark", required=True,
                        help="Specify the directory of the GenemarkS-T protein-peptides.csv file")
    parser.add_argument('-t', '--transdecoder', action="store", dest="transdecoder", required=True,
                        help="Specify the directory of the Transdecoder protein-peptides.csv file")
    parser.add_argument('-p', '--prefix', action="store", dest="prefix", default="sample",
                        help="Provide a prefix for the output csv files")
    args = parser.parse_args()

    try:
        os.makedirs("output/comparison_output")
    except FileExistsError:
        pass

    try:
        print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        genemark_data = csv_dataframe.extract_csv_data(args.genemark)
        transdecoder_data = csv_dataframe.extract_csv_data(args.transdecoder)

        find_distinct_peptides(transdecoder_data, genemark_data, args.prefix)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError:
        print(__doc__)
        print("Please provide valid files")
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
