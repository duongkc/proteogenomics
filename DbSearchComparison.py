#!/usr/bin/python3
import sys
import datetime
import re
import pandas

"""
Short script to compare the PEAKS psm search output to find overlap and distinction between the matched peptides.
"""


def clean_peptide_col(peptide_column):
    """Cleans up the peptide column by removing unnecessary information and returns the peptide"""
    sample_line = "K.ELC(+57.02)EQ(+.98)EC(+57.02)EWEEITITGSDGSTR.V"
    no_parentheses_pep = re.sub(r'\([^()]*\)', '', peptide_column)
    stripped_pep = no_parentheses_pep.replace('.', '')
    return stripped_pep


def extract_csv_data(input_file, db_type):
    """Parses the PEAKS protein-peptide csv file to extract the peptides and matching ORF accessions"""
    # And write them to a new temporary file, (should be split into sub methods tbh
    if db_type == "genemark":
        temp_output = "output/temp_gm_output.csv"
    else:
        temp_output = "output/temp_td_output.csv"

    # with open(temp_output, "w+") as temp_peptide_output:
    #     temp_peptide_output.write("Protein Accession,Peptide\n")
    #     with open(input_file, "r") as f:
    #         next(f)
    #         for line in f:
    #             columns = line.split(",")
    #             accession = columns[2]
    #             peptide_col = columns[3]
    #             peptide = clean_peptide_col(peptide_col)
    #             new_line = accession + "," + peptide + "\n"
    #             temp_peptide_output.write(new_line)
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    for i, row in csv_data.iterrows():
        raw_peptide = csv_data.at[i, 'Peptide']
        csv_data.at[i, 'Peptide'] = clean_peptide_col(raw_peptide)
    with open(temp_output, "w+") as temp_peptide_output:
        csv_data.to_csv(temp_peptide_output, sep=',', mode='w', index=False)



def find_distinct_peptides(transdecoder_file, genemark_file):
    """Filters the CSV files so only distinct peptides remain"""
    with open("output/distinct_td.csv", "w+") as distinct_transdecoder, \
            open("output/distinct_gm.csv", "w+") as distinct_genemark:
        transdecoder_data = pandas.read_csv(transdecoder_file, header='infer').drop_duplicates()
        genemark_data = pandas.read_csv(genemark_file, header='infer').drop_duplicates()
        td_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='left', indicator=True) \
            .query("_merge == 'left_only'")
        gm_merged = pandas.merge(transdecoder_data, genemark_data, on='Peptide', how='right', indicator=True) \
            .query("_merge == 'right_only'")
        td_merged[['Protein Accession_x', 'Peptide']] \
            .to_csv(distinct_transdecoder, sep=',', mode='w', index=False, header=['Protein Accession', 'Peptide'])
        gm_merged[['Protein Accession_y', 'Peptide']] \
            .to_csv(distinct_genemark, sep=',', mode='w', index=False, header=['Protein Accession', 'Peptide'])


extract_csv_data("data/propep_g.csv", "genemark")
# extract_csv_data("data/propep_t.csv", "transdecoder")
# find_distinct_peptides("output/temp_td_output.csv", "output/temp_gm_output.csv")


def main():
    print("Hello, World!")
