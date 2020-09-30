#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thurs Jun 4, 2020
Updated on

Vasu Kuram & Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

csi_to_gris_hgsc.py
    Reads the sample mapping and VCF information and 
    preps the metadata for uploading to GRIS

Quality Checks:
    *  which individuals have been returned in this batch
    *  each individual belongs to only one family
    *  Order IDs are consistent between pedigree, sample mapping, sample key and orders sent
    *  Finds and removes the duplicate and the control entries
    
"""
__author__ = 'Vasu Kuram & Susan Huse'
__date__ = 'June 4, 2020'
__version__ = '1.1'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input
import subprocess
import json
import fnmatch
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query, return_bsi_info
from generate_seqr_ped import generate_genrptlinks_ped


# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def query_bsi(ids, fields, query_field):
    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    bsi = bsi_query(curl_get, url_reports, session, fields, ids, query_field)

    return bsi


# batch_num = Batch Number (i.e. 23)
# unsequenced_fam = family members of patients received in current batch (before splitting) that do not have data received/processed yet
# masterkey = masterkey file for one half of the current batch
# received = all samples pulled straight from the sample key file delivered by HGSC
# moved_samples = list of sample IDs added to current batch that were moved from previous batches into current
def make_sample_tracking_file(batch_num, unsequenced_fam, masterkey, received, moved_samples):

    masterkey = masterkey[['Phenotips_Family_ID', 'Phenotips_ID', 'Exome_ID', 'DLM_LIS_Number', 'CRIS_Order#', 'Batch_Sent', 'Batch_Received']]

    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    batch_name = 'BATCH' + str(batch_num)
    batch_name_lower = 'Batch ' + str(batch_num)

    f = open('sample_tracking_summary_batch' + str(batch_num) + '.txt', 'w')
    f.write('******************************************\n*\n* Notes for samples released with Batch ' + str(batch_num) + '\n' +
            '*\n******************************************\n\n')

    fields = ['Phenotips ID', 'CRIS Order #', 'Batch Sent', 'Batch Received', 'CRIS Order Status', 'Active Status']
    sent = bsi_query(curl_get, url_reports, session, fields, [batch_name], 'Batch Sent')
    sent = sent[~sent['CRIS_Order_Status'].str.contains('cancel', case = False, na = False)]
    sent = sent[~sent['CRIS_Order_Status'].str.match('auto complete', case = False, na = False)]
    sent = sent[sent['Active status'].str.match('Active')]
    sent.drop_duplicates(keep = 'first', inplace = True)
    num_sent = len(sent['CRIS_Order#'].unique())

    f.write('Total number of samples sent to HGSC with Batch ' + str(batch_num) + ': ' + str(num_sent) + '\n')
    # f.write('** Samples sent together in ' + batch_name_lower + ' were split into two separate batches when received **\n\n')

    num_released_all = received.shape[0]
    f.write('Total number of new samples in latest HGSC release (before split): ' + str(num_released_all) + '\n\n')

    num_released = masterkey[masterkey['Batch_Received'].str.match(batch_name, na = False)].shape[0]
    f.write('Number of new samples released from HGSC in ' + batch_name_lower + ': ' + str(num_released) + '\n\n')

    num_in_masterkey = masterkey.shape[0]
    f.write('Total number of samples in the masterkey file: ' + str(num_in_masterkey) + "\n")
    f.write("Here's the breakdown: \n")
    f.write('---------------------------------------------------------------------------------------------\n\n')

    num_released_from_sent = len(intersection(sent['CRIS_Order#'].values.tolist(), masterkey['CRIS_Order#'].values.tolist()))
    f.write(str(num_released_from_sent) + ' sample(s) sent to HGSC in ' + batch_name_lower + ' were released with ' + batch_name_lower + '.\n\n')

    ##### Special Case: remove samples that are being added/moved from other batches into current batch from previous batches.
    # We'll make a separate section to outline those samples

    sent_in_other_and_released = masterkey[~masterkey['Batch_Sent'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[sent_in_other_and_released['Batch_Received'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[~sent_in_other_and_released['Phenotips_ID'].isin(moved_samples)]
    f.write(str(sent_in_other_and_released.shape[0]) + ' sample(s) sent to HGSC in other batches were released with ' + batch_name_lower + ':\n')
    f.write(sent_in_other_and_released.to_string(index = False) + '\n\n')

    # This section outlines samples moved to current batch:
    num_moved = len(moved_samples)
    if num_moved > 0:
        f.write(str(num_moved) + ' sample(s) released in previous batches were moved to ' + batch_name_lower + ':\n')
        f.write(masterkey[masterkey['Phenotips_ID'].isin(moved_samples)].to_string(index = False) + '\n\n')

    # all family members in the masterkey file (they're sequenced and have data on file)
    fam_members_released = masterkey[~masterkey['Batch_Received'].str.match(batch_name)]
    fam_members_released = fam_members_released[~fam_members_released['Phenotips_ID'].isin(moved_samples)]
    # fam_members_released = fam_members_released[['']]
    if fam_members_released.empty:
        f.write('0 sample(s) released in previous batches are family members of sample(s) released with ' + batch_name_lower + '.\n\n')
    else:
        f.write(str(fam_members_released.shape[0]) + ' sample(s) released in previous batches are family members of sample(s) released in ' + batch_name_lower + ':\n')
        f.write(fam_members_released.to_string(index = False) + '\n\n')

    f.write('---------------------------------------------------------------------------------------------\n\n')

    # get any unsequenced family members that are related to people in this half of the batch
    fam_members_unreleased = unsequenced_fam[unsequenced_fam['Phenotips_Family_ID'].isin(masterkey['Phenotips_Family_ID'])]
    if fam_members_unreleased.empty:
        f.write('There are no unreleased family members of sample(s) released in ' + batch_name_lower + ' on file.\n\n')
    else:
        f.write(str(fam_members_unreleased.shape[0]) + ' family members have been consented but not yet released.\n')
        fam_members_unreleased = fam_members_unreleased[['Phenotips_Family_ID', 'Phenotips_ID']]
        f.write(fam_members_unreleased.to_string(index = False))
        f.write('\n\n')

    sent_not_released = sent[~sent['CRIS_Order#'].isin(masterkey['CRIS_Order#'])]
    if sent_not_released.empty:
        f.write('All non-canceled samples sent to HGSC in ' + batch_name_lower + ' have been released in ' + batch_name_lower + '\n\n')
    else:
        f.write(str(sent_not_released.shape[0]) + ' sample(s) sent to HGSC in ' + batch_name_lower + ' were not released in ' + batch_name_lower + ':\n')
        f.write(sent_not_released.to_string(index=False) + '\n\n')

    f.close()


# Write shell script to symbolically link old bam files to new directory
# masterkey: one of the split masterkeys
def write_bam_link_script(masterkey, batch_num):
    batch_name = 'BATCH' + str(batch_num)
    previous_batches = masterkey[~masterkey['Batch_Received'].str.match(batch_name)]

    # Write script to link bams from previous batches (family/added samples) on LOCUS
    locus_dir = '/hpcdata/dir/CSI_DATA_PROCESSED'
    fname = 'link_previous_bams_' + 'batch' + str(batch_num) + '.sh'

    script = open(fname, "w")
    script.write("#!/bin/sh\nset -e\n\n")

    for i in previous_batches.index.values.tolist():
        pbatchdir = previous_batches.loc[i]['Batch_Received']
        pbatchdir = pbatchdir.replace('BATCH0', 'BATCH')
        phen_id = previous_batches.loc[i]['Phenotips_ID']

        bamlink = "ln -s {}/{}/BAM/{}.bam {}/{}/BAM/\n".format(locus_dir, pbatchdir, phen_id, locus_dir,
                                                               'BATCH' + str(batch_num))
        bailink = "ln -s {}/{}/BAM/{}.bam.bai {}/{}/BAM/\n".format(locus_dir, pbatchdir, phen_id, locus_dir,
                                                                   'BATCH' + str(batch_num))
        script.write(bamlink)
        script.write(bailink)
        script.write("\n")

    script.close()


def write_split_masterkeys(masterkey):

    first_half_name = 'BATCH' + str(batch)
    second_half_name = 'BATCH' + str(batch + 1)

    # sort df by date of enrollment
    masterkey['Date of Enrollment'] = masterkey['Date of Enrollment'].str.split(" ").str[0]
    masterkey['Date of Enrollment'] = pd.to_datetime(masterkey['Date of Enrollment'], format="%Y-%m-%d")
    masterkey = masterkey.sort_values(by=['Date of Enrollment'])
    masterkey.reset_index(drop=True, inplace=True)

    # divide the size in half to find midpoint
    midpoint = int(masterkey.shape[0] / 2)
    first_half = masterkey[:midpoint]
    second_half = masterkey[midpoint:]

    # find shared families
    shared_fams = intersection(first_half.Phenotips_Family_ID.values.tolist(),
                               second_half.Phenotips_Family_ID.values.tolist())

    # Remove shared families from first and second halves (we'll redistribute the families later)
    first_half = first_half[~first_half['Phenotips_Family_ID'].isin(shared_fams)]
    second_half = second_half[~second_half['Phenotips_Family_ID'].isin(shared_fams)]

    shared_fams = pd.DataFrame(shared_fams, columns=['Phenotips_Family_ID'])
    shared_fams['Min_Date'] = None
    shared_fams.drop_duplicates(keep='first', inplace=True)
    shared_fams.reset_index(drop=True, inplace=True)

    # get earliest date of a family member for each family
    for i in range(0, shared_fams.shape[0]):
        shared_fams.at[i, 'Min_Date'] = min(
            masterkey[masterkey['Phenotips_Family_ID'] == shared_fams.at[i, 'Phenotips_Family_ID']][
                'Date of Enrollment'])

    # put shared families in between the first and second half and re-split
    shared_fam_info = masterkey[masterkey['Phenotips_Family_ID'].isin(shared_fams['Phenotips_Family_ID'])]
    shared_fam_info = shared_fam_info.merge(shared_fams, how='inner', on='Phenotips_Family_ID')
    shared_fam_info = shared_fam_info.sort_values(by=['Min_Date', 'Phenotips_Family_ID'])
    shared_fam_info.reset_index(drop=True, inplace=True)
    # print(shared_fam_info[['Phenotips_Family_ID', 'Min_Date']])
    shared_fam_info.drop(columns=['Min_Date'], inplace=True)
    # print(shared_fam_info)

    masterkey_new = first_half.append(shared_fam_info, ignore_index=True)
    masterkey_new = masterkey_new.append(second_half, ignore_index=True)
    masterkey.reset_index(drop=True, inplace=True)

    # re-split
    first_half = masterkey_new[:midpoint]
    second_half = masterkey_new[midpoint:]

    # pull in any family members that might have been left
    # out from the split in the last family at the end of first_half
    if second_half.Phenotips_Family_ID.str.contains(first_half.at[midpoint - 1, 'Phenotips_Family_ID']).any():
        fam_id = first_half.at[midpoint - 1, 'Phenotips_Family_ID']
        fam_data = second_half[second_half.Phenotips_Family_ID == fam_id]
        second_half = second_half[~second_half.Phenotips_Family_ID.str.match(fam_id, na=False)]
        first_half = first_half.append(fam_data, ignore_index=True)

    first_half.reset_index(drop=True, inplace=True)
    second_half.reset_index(drop=True, inplace=True)

    first_half.drop(columns='Date of Enrollment', inplace=True)
    second_half.drop(columns='Date of Enrollment', inplace=True)

    second_half.Batch_Received.replace(first_half_name, second_half_name, inplace=True)

    # Test if anything is wrong
    shared_samples = intersection(first_half.Phenotips_Family_ID.values.tolist(), second_half.Phenotips_Family_ID.values.tolist())
    if len(shared_samples) == 0:
        print('\nSuccessfully split batches: ')
        print(first_half_name + ": " + str(first_half.shape[0]) + " samples")
        print(second_half_name + ": " + str(second_half.shape[0]) + " samples")
    else:
        print('Error occurred in batch splitting: Family members are still in both batches!')
        sys.exit()

    first_half.sort_values(by=['Batch_Received'], ascending = False, inplace = True)
    second_half.sort_values(by=['Batch_Received'], ascending = False, inplace = True)

    return first_half, second_half


# Reads in Sample Key file from HGSC
# Makes unsequenced_family df: all family members of current batch who have Batch_Received == blank
# Makes sequenced_family df: all family members of current batch who have Batch_Received filled in
def write_masterkey(sample_key_path):

    # read in sample key provided by Baylor and isolate the LIS number and exome ID (INDEX ID)
    received = pd.read_excel(sample_key_path)
    received = received[['INDEX ID', 'COLLABORATOR SAMPLE ID', 'FLOWCELL ID', 'LANE NUM']]
    received.rename(columns={'COLLABORATOR SAMPLE ID': 'DLM_LIS_Number'}, inplace=True)

    print('------------------------- Raw Sample Key -------------------------')
    print(str(received.shape[0]) + ' total samples received.')

    na12878 = received[received['DLM_LIS_Number'].str.contains('NA12878', na = False)]
    received = received[~received['DLM_LIS_Number'].str.contains('NA12878', na = False)]

    print('- ' + str(len(na12878['DLM_LIS_Number'])) + ' NA12878 control sample(s)')
    print('- ' + str(received.shape[0]) + ' actual samples\n\n')

    # make sample ID out of concatenation
    received['Sample_ID'] = received['FLOWCELL ID'].astype(str) + '-' + received['LANE NUM'].astype(str) + '-' + received['INDEX ID']
    received = received[['Sample_ID', 'DLM_LIS_Number']]

    # first get all family member IDs
    family_ids = query_bsi(received['DLM_LIS_Number'], ['Phenotips Family ID'], 'DLM LIS Number')
    family_ids = family_ids['Phenotips_Family_ID'].unique()
    # pd.DataFrame(family_ids).to_csv('fam_ids.csv')

    # Get IDs to merge from BSI
    bsi_fields = ['Phenotips ID', 'Batch Sent', 'Batch Received', 'DLM LIS Number', 'Exome ID',
                  'CRIS Order Status', 'CRIS Order #', 'Archive', 'Active Status', 'Phenotips Family ID', 'Date of Enrollment', 'Vendor']
    received_and_fam = query_bsi(family_ids, bsi_fields, 'Phenotips Family ID')
    # received_and_fam.to_csv('/Users/kuramvs/Documents/hgsc_scripts/test_files/id_dict.csv')

    # Clean up BSI data
    received_and_fam = received_and_fam[~received_and_fam['CRIS_Order_Status'].str.contains('cancel', case = False, na = False)]
    received_and_fam = received_and_fam[~received_and_fam['CRIS_Order_Status'].str.match('auto complete', case = False, na = False)]
    inactive = received_and_fam[~received_and_fam['Active status'].str.match('Active', case = False, na = False)]

    received_and_fam.drop_duplicates(keep = 'first', inplace = True)
    received_and_fam.replace('', np.NaN, inplace = True)
    received_and_fam.rename(columns={'Exome ID': 'Exome_ID', 'DLM LIS Number': 'DLM_LIS_Number'}, inplace=True)

    if inactive.shape[0] > 0:
        print('Warning! The following samples are Inactive, but were delivered with the latest CIDR release: ')
        print(inactive)
        print('')

    # print(received_and_fam.columns)
    # received_and_fam.to_csv('id_dict_cleaned.csv')

    # Merge additional BSI fields and rearrange columns
    received = received.merge(received_and_fam, on = 'DLM_LIS_Number', how = 'left')
    received = received[['Sample_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received', 'CRIS_Order#', 'Phenotips_Family_ID', 'Date of Enrollment', 'Vendor']]
    received.rename(columns = {'Sample_ID': 'Exome_ID'}, inplace = True)
    received['Batch_Received'] = "BATCH" + str(batch)

    received_and_fam = received_and_fam[received_and_fam['Active status'].str.match('Active', case = False, na = False)]
    received_and_fam.drop_duplicates(keep='first', inplace=True)

    # Get sequenced family members
    all_family = received_and_fam[~received_and_fam['Phenotips_ID'].isin(received['Phenotips_ID'])]

    sequenced_family = all_family[~all_family['Batch_Received'].isna()]
    sequenced_family = all_family[~all_family['Exome_ID'].isna()]
    sequenced_family = sequenced_family[['Exome_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received', 'Phenotips_Family_ID', 'Date of Enrollment', 'Vendor']]
    unsequenced_family = all_family[all_family['Batch_Received'].isna()]
    unsequenced_family = unsequenced_family[~unsequenced_family['Phenotips_ID'].isin(sequenced_family['Phenotips_ID'])]

    # add sequenced family to masterkey
    masterkey = received.append(sequenced_family, sort = False)

    print('----------------------- Un-split Masterkey -----------------------')
    print(str(masterkey.shape[0]) + ' total')
    print('- ' + str(received.shape[0]) + ' received')
    print('- ' + str(sequenced_family.shape[0]) + ' sequenced family\n\n')

    # Remove space from Batch Sent in BSI
    masterkey.Batch_Sent = masterkey.Batch_Sent.str.replace(' ', '')
    sequenced_family.Batch_Sent = sequenced_family.Batch_Sent.str.replace(' ', '')
    unsequenced_family.Batch_Sent = unsequenced_family.Batch_Sent.str.replace(' ', '')
    received.Batch_Sent = received.Batch_Sent.str.replace(' ', '')
    
    return masterkey, sequenced_family, unsequenced_family, received


# writes table containing Phenotips_ID, Exome_ID, Batch_Received, CRIS_Order#
# of all samples received from HGSC for Xi to upload to BSI (one file for each split batch)
def write_bsi_info(first_half, second_half):

    # Isolate only the new samples
    bsi_first = first_half[first_half['Batch_Received'].str.match('BATCH' + str(batch), na = False)]
    bsi_first = bsi_first[['Phenotips_ID', 'Exome_ID', 'Batch_Received', 'CRIS_Order#']]
    bsi_first['Vendor'] = 'HGSC'
    bsi_first.rename(columns = {'CRIS_Order#': 'CRIS_Order_ID'}, inplace = True)
    bsi_first.to_csv('Batch' + str(batch) + '_BaylorExoIDs.csv', index = False)

    if second_half is not None:
        bsi_second = second_half[second_half['Batch_Received'].str.match('BATCH' + str(batch + 1), na = False)]
        bsi_second = bsi_second[['Phenotips_ID', 'Exome_ID', 'Batch_Received', 'CRIS_Order#']]
        bsi_second['Vendor'] = 'HGSC'
        bsi_second.rename(columns = {'CRIS_Order#': 'CRIS_Order_ID'}, inplace = True)
        bsi_second.to_csv('Batch' + str(batch + 1) + '_BaylorExoIDs.csv', index = False)


def main():

    pd.set_option('mode.chained_assignment', None)
    global batch

    #
    # Usage statement
    #
    parseStr = 'Reads a list of available files, performs quality control on the data,\n\
    and outputs the files needed for GRIS.\n\n\
    Usage:\n\
        csi_to_gris_hgsc.py -b batch -d raw_data_dir -s sample_key_file\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--batch', required=True, type=int, help='Batch number (integer)')
    parser.add_argument('-d', '--dir', required=True, type=str, help='Directory containing raw BAMs and gVCFs from HGSC')
    parser.add_argument('-s', '--sample_key', required=True, type=str,  help='Path to HGSC Sample Key file')
    parser.add_argument('-u', '--unsplit', required=False, action='store_true', help='Will not split the batch into two')

    args = parser.parse_args()
    batch = args.batch
    dir = args.dir
    sample_key_path = args.sample_key
    unsplit = args.unsplit

    print('\nReading in sample key file and making masterkey file...')
    # Get masterkey, and information on sequenced/unsequenced family members of received patients
    masterkey, sequenced_family, unsequenced_family, received = write_masterkey(sample_key_path)
    # print(unsequenced_family.columns)

    if unsplit:
        batch_name = 'BATCH' + str(batch)
        print('Creating files for unsplit ' + batch_name)

        write_bsi_info(masterkey, None)
        make_sample_tracking_file(batch, unsequenced_family, masterkey, received, [])
        write_bam_link_script(masterkey, batch)

        masterkey[['Exome_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received']].to_csv('masterkey_' + batch_name + '.txt', index=False, sep='\t')

    else:
        # split the main masterkey in half
        # these will be used to subset the VCF into two batches
        first_half, second_half = write_split_masterkeys(masterkey)

        batch_name_first = 'BATCH' + str(batch)
        batch_name_second = 'BATCH' + str(batch + 1)

        # Check if all the numbers match up
        if first_half.shape[0] + second_half.shape[0] != masterkey.shape[0]:
            print('WARNING: not all samples from the main masterkey are in the split masterkeys')

        print('----------------------- Batch ' + str(batch) + ' Masterkey -----------------------')
        print(str(first_half.shape[0]) + ' total')
        print('- ' + str(first_half[first_half['Batch_Received'] == batch_name_first].shape[0]) + ' received')
        print('- ' + str(first_half[first_half['Batch_Received'] != batch_name_first].shape[0]) + ' sequenced family\n')

        print('----------------------- Batch ' + str(batch + 1) + ' Masterkey -----------------------')
        print(str(second_half.shape[0]) + ' total')
        print('- ' + str(second_half[second_half['Batch_Received'] == batch_name_second].shape[0]) + ' received')
        print('- ' + str(second_half[second_half['Batch_Received'] != batch_name_second].shape[0]) + ' sequenced family\n\n')

        # Make files to link patient IDs with Batch Received in BSI
        write_bsi_info(first_half, second_half)

        print("\nGenerating sample tracking files...")
        # Write sample tracking summary for both batches
        make_sample_tracking_file(batch, unsequenced_family, first_half, received, [])
        make_sample_tracking_file(batch + 1, unsequenced_family, second_half, received, [])

        # Make the link_previous_bams files
        write_bam_link_script(first_half, batch)
        write_bam_link_script(second_half, batch + 1)

        # Finally, write out all the masterkeys in desired format
        masterkey[['Exome_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received']].to_csv('masterkey.txt', index=False, sep='\t')
        first_half[['Exome_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received']].to_csv('masterkey_' + 'BATCH' + str(batch) + '.txt', index=False, sep='\t')
        second_half[['Exome_ID', 'Phenotips_ID', 'DLM_LIS_Number', 'Batch_Sent', 'Batch_Received']].to_csv('masterkey_' + 'BATCH' + str(batch + 1) + '.txt', index=False, sep='\t')


    print('\nAll done!')


if __name__ == '__main__':
    main()

