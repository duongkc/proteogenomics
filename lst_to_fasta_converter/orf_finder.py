#!/usr/bin/python3
"""short python module to extract ORFs found with GeneMarkS-T from the RNA transcripts assembled by Trinity.
"""

import argparse
import datetime
import os
import sys

nuc_dict = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}


def find_fasta_line(trans_file, accession):
    """Finds the right RNA transcript using the TRINITY gene accession found in the lst file"""
    with open(trans_file, "r") as inF:
        for line in inF:
            if accession in line:
                if line.startswith(">"):
                    sequence = next(inF)
                    return sequence
    return 0


def extract_gene_id(line):
    """Extracts gene accession of the target transcript"""
    gene_id = line.split()[3]

    return gene_id


def find_transcript_position(info_line):
    """Finds strand and position within the target transcripts of the final selected ORFs"""
    columns = info_line.split()
    strand = columns[1]
    start = int(columns[2].strip("<")) - 1
    stop = int(columns[3].strip(">"))
    return strand, start, stop


def reverse_complement(sequence):
    """Converts sequence to reverse complement"""
    bases = list(sequence)
    bases = reversed([nuc_dict.get(base, base) for base in bases])
    new_sequence = ''.join(bases)
    return new_sequence


def write_orfs(gene_id, sequence, start, stop, strand):
    """Writes ORFs to new file in current directory"""
    with open("output/Trinity.fasta.genemark.cds", "a+") as new_file:
        new_file.write(">" + gene_id + ":" + str(start+1) + "-" + str(stop) + "(" + strand + ")\n")
        new_file.write(sequence + "\n")


def insert_newlines(seq, every=60):
    """adds \n after every 60 nucs in sequence"""
    lines = []
    for i in range(0, len(seq), every):
        lines.append(seq[i:i + every])
    return '\n'.join(lines)


def parse_genemark(gm_file, trans_file):
    """Parses through genemark file to extract info for each transcript"""
    with open(gm_file, "r") as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if 'FASTA' in line:
                gene_id = extract_gene_id(line)
                sequence = find_fasta_line(trans_file, gene_id)
                try:
                    info_line = lines[index + 4]
                    if info_line.strip():
                        strand, start_pos, stop_pos = find_transcript_position(info_line)
                        sequence = sequence[start_pos:stop_pos]
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        sequence = insert_newlines(sequence)
                        write_orfs(gene_id, sequence, start_pos, stop_pos, strand)

                except IndexError:
                    pass


def main(argv):
    print(' '.join(argv))

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-g', '--genemark', action="store", dest="genemark", required=True,
                        help="Specify the directory of the GenemarkS-T protein-peptides.csv file")
    parser.add_argument('-t', '--transcript', action="store", dest="transcript", required=True,
                        help="Specify the directory of the Trinity RNA transcript fasta file")
    args = parser.parse_args()

    try:
        os.makedirs("output")
    except FileExistsError:
        pass

    if os.path.exists("output/Trinity.fasta.genemark.cds"):
        os.remove("output/Trinity.fasta.genemark.cds")

    try:
        print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
        parse_genemark(args.genemark, args.transcript)
        print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    except FileNotFoundError:
        print(__doc__)
        print("Please provide valid files")


if __name__ == '__main__':
    main(sys.argv)
