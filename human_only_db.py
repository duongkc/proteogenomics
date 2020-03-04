#!/usr/bin/python3
"""
Filters out any non human entries from a fasta database
To be used to create a swissprot database with only Homo sapiens sequences
"""
import argparse
import datetime
import os
import sys

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
    try:
        print("Started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        search_db(args.database)
        print("Finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError as e:
        print(__doc__)
        print("Please provide valid files:")
        print(e)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)