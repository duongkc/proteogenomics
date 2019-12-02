#!/usr/bin/python3

import sys
import datetime
import pandas
import getopt
import os

def extract_comparison_csv_data():
    """Reads the comparison csv data as dataframe"""
    input_file = "comparison_output/test1_distinct_gm.csv"
    csv_data = pandas.read_csv(input_file, header='infer', delimiter=',')
    