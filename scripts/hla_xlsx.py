#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:09:26 2019


@author: Vasu Kuram
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

ncbr_bsi and ncbr_huse by Susan Huse

hla_table.py
    Generates HLA allele table containing HLA alleles for every subject in the batch.
    Switches from Phenotips ID to CIDR Exome ID and creates CSV files (one with Phenotips and one with CIDR)
    
Assumes that this script is run inside the directory for each Batch i.e. in BATCH15 folder
    
####################################
The following lines are to be run before executing this script:
    module load anaconda3/5.3.0
    module load curl/7.60.0-GCCcore-7.3.0

"""

import sys
import os
import re
import datetime
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import datetime
import time
import json
import xlsxwriter

# Simulates Excel autofit. Gets width for each column.
def get_col_widths(dataframe):
    widths = []

    for column in dataframe:
        strings = dataframe[column].values.tolist()
        strings.append(column)
        widths.append(len(max(strings, key=len)) + 1)

    return widths


# Removes locus name from allele type
def get_hla_type(hla_types):
    types = []

    for i, row in hla_types.iterrows():
        to_replace = row.Locus + '*'
        types.append(row.Allele.replace(to_replace, ''))

    return types


# Constructs and outputs an HLA Table for each family in the batch
def make_table(subject_id, input, output):
    path = '/hpcdata/dir/CIDR_DATA_RENAMED/'
    colnames = ['Subject_ID']

    hla_path = input
    hla_tab = pd.DataFrame()

    try:  # ID may be missing so test it
        df = pd.read_csv(hla_path, delimiter="\t", dtype=str)
    except:  # Writes out to error log and skips to next iteration
        print(hla_path)
        print('Error: Cannot open input file ' + input + '\n')
        return
    else:  # if no error occurs, append the row

        # concatenate/append IDs and hla types for allele 1
        allele1 = get_hla_type(df[df['Chromosome'] == '1'])
        row = pd.DataFrame(np.concatenate(([subject_id], allele1), axis=0)).T
        hla_tab = hla_tab.append(row, ignore_index=True)

        # concatenate/append IDs and hla types for allele 2
        allele2 = get_hla_type(df[df['Chromosome'] == '2'])
        row = pd.DataFrame(np.concatenate(([''], allele2), axis=0)).T
        hla_tab = hla_tab.append(row, ignore_index=True)

        if len(colnames) == 1:  # if no column names, grab from first input file
            loci = df['Locus'].unique()
            allele_cols = ["{}*".format(locus) for locus in loci]
            colnames = colnames + allele_cols
    #                     print(colnames)

        hla_tab.columns = colnames
        fname = output

        # Get the xlsxwriter workbook and worksheet objects.
        workbook = xlsxwriter.Workbook(fname)
        worksheet = workbook.add_worksheet(name='HLA')

        # Add a header format.
        header_format = workbook.add_format({
            'bold': True,
            'valign': 'top',
            'fg_color': '#bfbfbf',
            'border': 0,
            'bottom': 1})

        row_format1 = workbook.add_format({
            'fg_color': '#ffffff',
            'border': 0})
        row_format2 = workbook.add_format({
            'fg_color': '#ffffff',
            'border': 0,
            'bottom': 1})

        # Write the column headers with the defined format.
        for col_num, value in enumerate(hla_tab.columns.values):
            worksheet.write(0, col_num, value, header_format)

        for col in range(0, hla_tab.shape[1]):
            counter = 1
            for row in range(0, hla_tab.shape[0]):
                if counter == 1:
                    worksheet.write(row + 1, col, hla_tab.iloc[row, col], row_format1)
                    counter = counter + 1
                elif counter == 2:
                    worksheet.write(row + 1, col, hla_tab.iloc[row, col], row_format2)
                    counter = 1

        for i, width in enumerate(get_col_widths(hla_tab)):
            worksheet.set_column(i, i, width)

        worksheet.write(hla_tab.shape[0] + 2, 0, "NOTE:")
        worksheet.write(hla_tab.shape[0] + 3, 0,
                        "HLA types are not phased across loci, data rows do not represent results from a particular chromosome.")
#        worksheet.write(hla_tab.shape[0] + 4, 0, "Accuracy has not yet been evaluated for the 2nd or 3rd fields.")
        worksheet.write(hla_tab.shape[0] + 4, 0,
                        "The HLA-LA author assessed accuracy at G group resolution (fields 1, 2 and 3) for a subset of the 1000 Genomes data set and showed a 100% agreement with clinical HLA typing.")
        worksheet.write(hla_tab.shape[0] + 5, 0,
                        "Treat calls for HLA-DRB3/4 with caution. The model does not try to estimate copy number for these genes, so you'll end up with calls for genes that are not present")
        worksheet.write(hla_tab.shape[0] + 6, 0,
                        "(in interpreting calls for the DRB orthologs, it helps to take into account HLA-DRB1 genotype, which, generally speaking and via linkage, will tell you how many HLA-DRB3/4/5 copies to expect).")

        workbook.close()
        # Close the Pandas Excel writer and output the Excel file.


def main():
    # Usage statement and options
    parseStr = 'Creates .xlsx file using an R1 bestguess HLA .txt file \n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, help = 'Input R1 bestguess file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output .xlsx file name')
    parser.add_argument('-s', '--subject_id', required=True, type=str, help='Subject ID')

    args = parser.parse_args()
    input = args.input
    output = args.output
    subject_id = args.subject_id

    print('Generating tables...')
    make_table(subject_id, input, output)
    print('All done!')


if __name__ == '__main__':
    main()
