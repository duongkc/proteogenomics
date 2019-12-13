#!/usr/bin/python3
"""
Module that converts cds sequences to pep.
"""

import datetime


def translate(sequence):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein = ""
    seq = sequence.rstrip()
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


def insert_newlines(seq, every=60):
    """adds \n after every 60 nucs in sequence"""
    lines = []
    for i in range(0, len(seq), every):
        lines.append(seq[i:i + every])
    return '\n'.join(lines)


def parse_cds_fasta(file):
    """Extracts cds data to convert and write to pep fasta file"""
    with open("../output/Trinity.fasta.genemark.pep", "w+") as new_file:

        with open(file, "r") as f:
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    if sequence:
                        pep = translate(sequence)
                        fasta_pep = insert_newlines(pep)
                        new_file.write(fasta_pep + "\n")
                    sequence = ""
                    header = line
                    new_file.write(header)
                else:
                    sequence += line.rstrip()
            pep = translate(sequence)
            fasta_pep = insert_newlines(pep)
            new_file.write(fasta_pep + "\n")


def main():

    print("started at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))
    parse_cds_fasta("../output/Trinity.fasta.genemark.cds")
    print("finished at: " + datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))


if __name__ == '__main__':
    main()
