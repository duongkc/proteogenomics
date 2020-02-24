#!/usr/bin/python3
"""
Filters out any non human entries from a fasta database
To be used to create a swissprot database with only Homo sapiens sequences
"""
import argparse
import datetime
import multiprocessing as mp
import os
import sys

import numpy as np
from Bio import SeqIO


def search_db(db_file):
    filename = os.path.splitext(db_file)[0]
    output = "{}.human.fasta".format(filename)
    with open(output, "w") as out_file:
        with open(db_file, "r") as database:
            for record in SeqIO.parse(database, "fasta"):
                if "Homo sapiens" in record.description:
                    SeqIO.write(record, out_file, "fasta")


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', '--database', action='store', dest="database", required=True,
                        help="Specify the directory of the protein database file")
    args = parser.parse_args()

    search_db(args.database)


if __name__ == '__main__':
    main(sys.argv)