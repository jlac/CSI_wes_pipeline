#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:43:58 2018
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

csi_to_gris.py
    Reads the sample mapping and VCF information and 
    preps the metadata for uploading to GRIS

v 1.0 - initial code version.
  1.1 - moved all queries to automated BSI query
        updated to move 00 as mother phenotips to 0 for unknown
        renamed README to sample_tracking_summary
  1.2 - added line to remove CRIS Order Status = Auto Complete from family members query results

Quality Checks:
    *  which individuals have been returned in this batch
    *  each individual belongs to only one family
    *  Order IDs are consistent between pedigree, sample mapping, sample key and orders sent
    *  Finds and removes the duplicate and the control entries
    
"""
__author__ = 'Susan Huse'
__date__ = 'September 5, 2018'
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

####################################
# 
# Functions for processing data
#
####################################
# Test file exists or error out
def test_file(f):
    if not os.path.isfile(f):
        err_out("Error: unable to locate input file:  " + f, log)
    else:
        return(True)

# Create config dictionary of input files
def create_config(dir_info):
    #print(dir_info)
    filenames = os.listdir(dir_info)
    #print(filenames)

    config = dict()
    # Subject Sample Mapping file
    f = fnmatch.filter(os.listdir(dir_info), "Holland*SubjectSampleMappingFile.csv")
    if len(f) < 1:
        err_out('Error: unable to locate input Subject Sample Mapping file: "*SubjectSampleMappingFile.csv"\n', log)
    else: 
        config['mapping'] = os.path.join(dir_info, f[0])

    # Master Sample Key file
    f = fnmatch.filter(os.listdir(dir_info), "Holland*MasterSampleKey*.csv")
    if len(f) < 1:
        err_out('Error: unable to locate input Master Sample Key file: "*MasterSampleKey.csv"\n', log)
    else:
        config['samplekey'] = os.path.join(dir_info, f[0])

    # Pedigree file
    f = fnmatch.filter(os.listdir(dir_info), "Holland*Pedigree*.csv")
    if len(f) < 1:
        err_out('Error: unable to locate input Pedigree file: "*Pedigree*.csv"\n', log)
    else:
        config['pedigree'] = os.path.join(dir_info, f[0])

    # test for the three files, will error out if any are missing
    foundfilestext = ""
    for i in ['mapping', 'samplekey', 'pedigree']:
        test_file(config[i])
        foundfilestext += "\t" + i + ": " + config[i] + "\n"

    foundfilestext = "Successfully located {} input files:\n{}\n".format(str(len(config)), foundfilestext)

    # add to the config variables
    config['rootdir'] = '/hpcdata/dir/CIDR_DATA_RENAMED'
    config['family_errors'] = 'family_errors.txt'
    
    # output file names
    config['tracking'] = 'sample_tracking_summary_batch' + str(batch) + '.txt'
    config['masterkey'] = 'masterkey_batch' + str(batch) + '.txt'
    config['newpedigree'] = 'seqr_ped_batch' + str(batch) + '.txt'
    config['batchinfo'] = 'genrptlinks_batch' + str(batch) + '.txt'
    config['linkscript_fname'] = 'link_bams_batch' + str(batch) + '.sh'

    return(config, foundfilestext)

# Compare indices of two data series, return missing values
def compare_keys(x,y):
    # join outer, so missing values with be null values in the x or y column

    missing_in_x = list(set(y.index.values) - set(x.index.values))
    missing_in_y = list(set(x.index.values) - set(y.index.values))

    return(missing_in_x, missing_in_y)


# Python program to illustrate the intersection
# of two lists in most simple way
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


######################################
#   
# Functions for importing data
#
######################################

# Find the duplicate ID from the pedigree file
def find_the_duplicate(df, fformat):
    if fformat == "ped":
        # find the duplicate using the comments column with the word "duplicate"
        dupestr = df[df['Investigator Column 1'].str.contains('Duplicate', case=False, na=False)]['Subject_ID']
        # find the duplicate because it is the only subject ID > 9 characters
        dupelen = df.loc[df['Subject_ID'].str.len() > 9]['Subject_ID']

        # removes any collaborator data
        dupelen = dupelen[dupelen.str.startswith('002')]

        # If the two methods aren't the same answer, error out
        
        # send_update("Length: " + str(dupestr))
#         send_update("Length: " + str(dupelen))
        if not dupestr.equals(dupelen):
            err_out("Warning: unable to ascertain duplicate sample identifier.\n" + \
                    "{} was identified in the pedigree file as a duplicate in 'Investigator Column 1',\n" + \
                    "{} has Subject_ID length > 9 characters.\n" + \
                    "Exiting".format(dupestr, dupelen), log)

    elif fformat == "manifest":
        # find the duplicate because it is the only subject ID > 9 characters
        dupestr = df.loc[(df['Subject_ID'].str.len() > 9) & (df['Subject_ID'] != "CONTROL_ID")]['Subject_ID']

    if dupestr.size > 1:
        send_update("Warning: found more than one duplicate sample identifier: {}.".format(dupestr.tolist()), log)
    elif dupestr.size < 1:
        send_update("Warning: no duplicate sample identifier was found in the {} file.".format(fformat), log)
        return(None)
    else:
        send_update("Duplicate sample identifier: {}".format(dupestr.tolist()), log)

    return(dupestr.tolist())

# Find the control ID from the 
def find_the_control(x):
    
    # search the indices for anything with NA or with HG
    regexNA = re.compile('^NA')
    regexHG = re.compile('^HG')
    ctrl = [i for i in x.index.tolist() if regexNA.match(i)]
    ctrlHG = [i for i in x.index.tolist() if regexHG.match(i)]
    ctrl.extend(ctrlHG)

    send_update("Control sample identifier: {}".format(set(ctrl)), log)

    if len(ctrl) != 1:
        error_text = "\nError: found {} control sample identifier(s), expected only one.\n{}".format(str(len(ctrl))    ,ctrl) 
        error_text = error_text + "\nPlease enter 'y' to continue or 'q' to quit.\n"
        pause_for_input(error_text, 'y', 'q', log)
        ctrl = None

    else:
        ctrl = ctrl[0]

    return(ctrl)


# Import pedigree information from Excel file
def import_pedigree(f):
    test_file(f)
    
    # import the data
    # but some are csv and some are excel! 
    if bool(re.search('.csv$', f)):
        peds = pd.read_csv(f, dtype={'Family': object}, encoding='utf-8')
    elif bool(re.search('.xlsx*$', f)):
        peds = pd.read_excel(f)
    else:
        err_out("Error: Confused by pedigree file name {}, expecting *.csv or *.xlsx.\nExiting.".format(f), log)

    peds = peds[peds['Subject_ID'] != ""]
    
    ## Each set has a sequencing duplicate that should be removed
    theDupe = find_the_duplicate(peds, "ped")

    # create a series of batch information with Subject as the index
    peds = peds[['Subject_ID', 'Investigator Column 3']]
    peds = peds[peds.Subject_ID.notnull()].set_index('Subject_ID')
    peds = peds['Investigator Column 3'].str.replace("_.*$", "", regex=True)

    if theDupe != None:
        peds = peds.drop(theDupe, errors='ignore')
        
#     peds.to_csv('peds.csv')

    return(peds, theDupe)

# Import sample mapping information from csv file, export collaborator masterkey
def import_samplemapping(f, findControl):
    test_file(f)
    collaborators = None

    # import and create series with index
    mapping = pd.read_csv(f, header=0, encoding='utf-8')
#    collaborators = mapping[~mapping['SUBJECT_ID'].str.startswith('002')]
#    collaborators.to_csv('collaborator_masterkey_' + batch_name + '.txt', sep='\t', header=True, index=False)
    
    mapping.columns.values[0] = "Subject_ID"
    mapping.columns.values[1] = "Exome_ID"

    # Export collaborator IDs that don't match usual format to another masterkey file
    if not pedonly:
        collaborators = mapping[~mapping['Subject_ID'].str.startswith("002")]
        if ~collaborators.empty:
            collaborators = collaborators[~collaborators['Subject_ID'].str.contains("NA12878")]
            if ~collaborators.empty:
                collaborators.drop(['SAMPLE_SOURCE', 'SOURCE_SAMPLE_ID'], axis = 1, inplace = True)
                collaborators.to_csv('collaborator_masterkey_batch' + str(batch) + '.txt', sep='\t', header=True, index=False)
    
    mapping = mapping.set_index('Subject_ID')['Exome_ID']

    # find the Control, and then remove it and the dupe
    # Will read old ones too and you don't want to change the control
    if(findControl):
        ctrl = find_the_control(mapping)
        if ctrl is not None:
            mapping = mapping.drop([ctrl])
    else:
        ctrl=None

    return(mapping, ctrl, collaborators)

# Import sample key information from csv file
def import_samplekey(f, ctrl):
    test_file(f)

    # import and convert to series
    samplekey = pd.read_csv(f, header=0, encoding='utf-8')
    samplekey = samplekey.set_index('Subject_ID')
    samplekey = samplekey['LIMS_SampleId']

    if ctrl != None:
        samplekey = samplekey.drop([ctrl])
    #samplekey = samplekey.drop([dupe, ctrl])

    return(samplekey)

######################################
#   
# Functions for writing output files
#
######################################

# Create Sample Tracking file (old README file)
def create_tracking(sentDF, receivedDF, familyDF, filename):
    # Header
    stars = "******************************************"
    trackingtxt = "\n".join([stars, "*", "* Notes for samples released with " + batch_label, "*", stars, "",""])

    # Sent
    trackingtxt += "Total number of samples sent to CIDR with {}: {}\n\n".format(batch_label, str(sentDF.shape[0]))

    # Received
    trackingtxt += "Total number of samples released from CIDR with {}: {}\n\n".format(batch_label, str(receivedDF.shape[0]))

    # Sent and Received
    received_current = receivedDF.loc[receivedDF['Batch_Sent'] == batch_name].shape[0]
    trackingtxt += "{} sample(s) sent to CIDR in {} were released with {}.\n\n".format(str(received_current), batch_label, batch_label)

    # Sent not Received, Received not from ealier batches
    sentIDs = set(sentDF['CRIS_Order#'])
    receivedIDs = set(receivedDF['CRIS_Order#'])
                                 
    # VASU: this just takes out sent from received 
    # update this to any sent IDs - receivedIDs in current batch or earlier batches?
    
    no_blanks = sentDF[sentDF['Batch_Received'] != '']
    sent_and_received = set(no_blanks[no_blanks['Batch_Received'].str.slice(5).astype(int) <= batch]['CRIS_Order#'])    # Samples that were sent in CURRENT batch but received in EARLIER batch 

    sent_not_received = list(sentIDs - receivedIDs - sent_and_received)
#    print(sent_not_received)
    # VASU: this should also be updated?
    received_from_other = list(receivedIDs - sentIDs)
    
    # Print Sent not Received
    sent_not_received_df = sentDF[sentDF['CRIS_Order#'].isin(sent_not_received)]
    if no_blanks.empty:
        trackingtxt += "All non-canceled samples sent to CIDR in {} have been released in {}:\n\n".format(batch_label, batch_label)
    else:
        trackingtxt += "{} sample(s) sent to CIDR in {} have not been released in {}:\n{}\n\n".format(str(len(sent_not_received)), batch_label, batch_label, sent_not_received_df.to_string(index=False))

    # Print Received from other batches (not sent from current batch)
    # VASU: this should be fixed if "received_from_earlier" above is fixed
    received_not_sent_df = receivedDF[receivedDF['CRIS_Order#'].isin(received_from_other)]  #contains samples that were received in current batch but not sent out from this batch

#    print(received_not_sent_df.shape[0])
#    print(received_from_prev_df.shape[0] + received_from_future_df.shape[0])
    
    trackingtxt += "{} sample(s) sent to CIDR in other batches were released with {}:\n{}\n\n".format(str(received_not_sent_df.shape[0]), batch_label, received_not_sent_df.to_string(index=False))

    # Find family members, previously released (seq) and not yet sequenced (unseq)
    # remove anything that was received with this batch, set as sequenced
    familyDF['Exome_ID'].replace('', np.nan, inplace=True)
    
    #both_family_members contains family members that are neither received in current batch nor received in other previous/future batches
    both_family_members = familyDF.loc[~familyDF['CRIS_Order#'].isin(receivedIDs)]  #filter out family members that were received in current batch
    both_family_members = both_family_members.loc[~both_family_members['CRIS_Order#'].isin(received_from_other)]    #filter out family members that were received in previous/future batch

    # Remove from sequenced anything without a batch received or anything that was canceled
    seq_family_members = both_family_members.loc[~both_family_members['CRIS_Order_Status'].str.contains('Cancel', case=False)]
    ### SUE: double check that this doesn't lose important records elsewhere
    seq_family_members = seq_family_members.loc[~seq_family_members['CRIS_Order_Status'].str.contains('Auto Complete', case=False)]
    seq_family_members = seq_family_members.loc[~seq_family_members['Exome_ID'].str.contains('Cancel', case=False, na=True)]
    seq_family_members = seq_family_members.loc[familyDF['Batch_Received'] != ""]
    seq_family_members = seq_family_members.loc[familyDF['Batch_Received'].str.contains("BATCH")]
    seq_family_members = seq_family_members[seq_family_members['Batch_Received'].str.slice(5).astype(int) < batch]

    # Remove from sequenced anything without a CIDR Exome ID, shouldn't be any but just in case
    #seq_family_members.dropna(subset=['Exome_ID'], inplace=True)
    seq_family_members = seq_family_members.loc[seq_family_members.Exome_ID.notnull()]

    # Subset columns for both dataframes
    family_columns = ['CRIS_Order#', 'Phenotips_ID', 'Phenotips_Family_ID', 'Batch_Sent', 'Batch_Received']
    #seq_family_members = seq_family_members[family_columns]

    # Set unsequenced as all of the family members that are not in the sequenced set
    # then drop duplicates rows (canceled orders, etc)
    unseq_family_members = both_family_members.loc[~both_family_members.Phenotips_ID.isin(seq_family_members.Phenotips_ID)]
    unseq_family_members = unseq_family_members.loc[~unseq_family_members.Phenotips_ID.isin(receivedDF['Phenotips_ID'])]
    ped_fields = ['Phenotips_Family_ID', 
                   'Phenotips_ID', 
                   'Father_Phenotips_ID', 
                   'Mother_Phenotips_ID', 
                   'Gender', 
                   'Affected']
    unseq_family_members = unseq_family_members[ped_fields].drop_duplicates()

    # VASU: probably the worst offender:
    if seq_family_members.empty:
        trackingtxt += "No sample(s) released in previous batches are family members of sample(s) released with {}.\n\n".format(batch_label)
    else:
        trackingtxt += "{} sample(s) released in previous batches are family members of sample(s) released with {}.\n{}\n\n".format(str(seq_family_members.shape[0]), batch_label, seq_family_members[family_columns].to_string(index=False))
    
    trackingtxt += "{} family members have been consented but not yet released.\n{}\n\n".format(str(unseq_family_members.shape[0]), unseq_family_members.to_string(index=False))

    if not pedonly:
        f = open(filename,'w')
        send_update(trackingtxt, f, True)

    # familyDF.to_csv('familyDF.csv', index = False)
    # seq_family_members.to_csv('seq_family_members.csv', index = False)
    
    return(familyDF, seq_family_members, unseq_family_members)

# batch_num = Batch Number (i.e. 23)
# fams = all samples in current batch + all family members (basically pedigree with more info columns)
# masterkey = masterkey file for current batch
# moved_samples = list of sample IDs added to current batch that were moved from previous batches into current
def make_sample_tracking_file(batch_num, fams, masterkey, sample_key, moved_samples):

    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    batch_name = 'BATCH' + str(batch_num)
    batch_name_lower = 'Batch ' + str(batch_num)

    f = open('sample_tracking_summary_batch' + str(batch_num) + '.txt', 'w')
    f.write('******************************************\n*\n* Notes for samples released with Batch ' + str(batch_num) + '\n' +
            '*\n******************************************\n\n')

    fields = ['Phenotips ID', 'CRIS Order #', 'Batch Sent', 'Batch Received', 'CRIS Order Status', 'Active Status']
    sent = bsi_query(curl_get, url_reports, session, fields, ['BATCH23'], 'Batch Sent')
    sent = sent[~sent['CRIS_Order_Status'].str.contains('cancel', case = False, na = False)]
    sent = sent[~sent['CRIS_Order_Status'].str.match('auto complete', case = False, na = False)]
    sent = sent[sent['Active status'].str.match('Active')]
    sent.drop_duplicates(keep = 'first', inplace = True)
    num_sent = len(sent['CRIS_Order#'].unique())

    f.write('Total number of samples sent to CIDR with Batch ' + str(batch_num) + ': ' + str(num_sent) + '\n')
    # f.write('** Samples sent together in ' + batch_name_lower + ' were split into two separate batches when received **\n\n')

    sample_key = pd.read_csv(config['samplekey'])
    sample_key = sample_key[sample_key['Subject_ID'].str.len() < 10]
    sample_key = sample_key[~sample_key['Subject_ID'].str.match('NA12878')]
    num_released_all = sample_key.shape[0]
    f.write('Total number of new samples in latest CIDR release: ' + str(num_released_all) + '\n\n')

    num_released = masterkey[masterkey['Batch_Received'].str.match(batch_name, na = False)].shape[0]
    f.write('Number of new samples released from CIDR in ' + batch_name_lower + ': ' + str(num_released) + '\n\n')

    num_in_masterkey = masterkey.shape[0]
    f.write('Total number of samples in the masterkey file: ' + str(num_in_masterkey) + "\n")
    f.write("Here's the breakdown: \n")
    f.write('---------------------------------------------------------------------------------------------\n\n')

    num_released_from_sent = len(intersection(sent['CRIS_Order#'].values.tolist(), masterkey['CRIS_Order#'].values.tolist()))
    f.write(str(num_released_from_sent) + ' sample(s) sent to CIDR in ' + batch_name_lower + ' were released with ' + batch_name_lower + '.\n\n')

    ##### Special Case: remove samples that are being added/moved from other batches into current batch from previous batches.
    # We'll make a separate section to outline those samples

    sent_in_other_and_released = masterkey[~masterkey['Batch_Sent'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[sent_in_other_and_released['Batch_Received'].str.match(batch_name, na = False)]
    sent_in_other_and_released = sent_in_other_and_released[~sent_in_other_and_released['Phenotips_ID'].isin(moved_samples)]
    f.write(str(sent_in_other_and_released.shape[0]) + ' sample(s) sent to CIDR in other batches were released with ' + batch_name_lower + ':\n')
    f.write(sent_in_other_and_released.to_string(index = False) + '\n\n')

    # This section outlines samples moved to current batch:
    num_moved = len(moved_samples)
    if num_moved > 0:
        f.write(str(num_moved) + ' sample(s) released in previous batches were moved to ' + batch_name_lower + ':\n')
        f.write(masterkey[masterkey['Phenotips_ID'].isin(moved_samples)].to_string(index = False) + '\n\n')

    fam_members_released = masterkey[~masterkey['Batch_Received'].str.match(batch_name, na = False)]
    fam_members_released = fam_members_released[~fam_members_released['Phenotips_ID'].isin(moved_samples)]
    if fam_members_released.empty:
        f.write('No sample(s) released in previous batches are family members of sample(s) released with ' + batch_name_lower + '.\n\n')
    else:
        f.write(str(fam_members_released.shape[0]) + ' sample(s) released in previous batches are family members of sample(s) released in ' + batch_name_lower + ':\n')
        f.write(fam_members_released.to_string(index = False) + '\n\n')

    f.write('---------------------------------------------------------------------------------------------\n\n')

    fam_members_unreleased = fams[~fams['Phenotips_ID'].isin(masterkey['Phenotips_ID'])]
    if fam_members_unreleased.empty:
        f.write('There are no unreleased family members of sample(s) released in ' + batch_name_lower + ' on file.\n\n')
    else:
        f.write(str(fam_members_unreleased.shape[0]) + ' family members have been consented but not yet released.\n')
        fam_members_unreleased = fam_members_unreleased[['Phenotips_Family_ID', 'Phenotips_ID', 'Father_Phenotips_ID', 'Mother_Phenotips_ID', 'Gender', 'Affected']]
        f.write(fam_members_unreleased.to_string(index = False))
        f.write('\n\n')

    sent_not_released = sent[~sent['CRIS_Order#'].isin(masterkey['CRIS_Order#'])]
    if sent_not_released.empty:
        f.write('All non-canceled samples sent to CIDR in ' + batch_name_lower + ' have been released in ' + batch_name_lower + '\n\n')
    else:
        f.write(str(sent_not_released.shape[0]) + ' sample(s) sent to CIDR in ' + batch_name_lower + ' were not released in ' + batch_name_lower + ':\n')
        f.write(sent_not_released.to_string(index=False) + '\n\n')

    f.close()


# Output the three additional files
def write_files(filenames, familyDF, receivedIDs, seqFamDF, sample_key, batch):
    #filenames 0: masterkey, 1: batchinfo, 2: seqrped
    receivedDF = familyDF.loc[familyDF['CRIS_Order#'].isin(receivedIDs)]

    #seqFamDF = familyDF[familyDF['Phenotips_ID'].isin(seqIDs)]

    receivedDF.Batch_Received = receivedDF.Batch_Received.replace('', np.nan)

    if receivedDF.Batch_Received.isna().any():
        receivedDF.Batch_Received = "BATCH" + str(batch)

    seqDF = pd.concat([receivedDF, seqFamDF])
    master_fields = ['CRIS_Order#',
                     'Phenotips_ID',
                     'Exome_ID',
                     'Batch_Sent',
                     'Batch_Received']
    masterkey = seqDF[master_fields]

    # Write masterkey: Contains sequenced family members + those received in current batch
    if not pedonly:
        masterkey.to_csv(filenames[0], sep='\t', header=True, index=False)

    # replace old pedigree generation method with BSI only method:
    ped, genrptlinks, all_fams = generate_genrptlinks_ped(batch)    #all fams contains same IDs as ped but with all BSI fields
    ped.to_csv(filenames[2], sep='\t', header=False, index=False)

    if not pedonly:
        genrptlinks.to_csv(filenames[1], sep='\t', header=False, index=False)
        make_sample_tracking_file(batch, all_fams, masterkey, sample_key, [])


# Write shell script to symbolically link old bam files to know directory
def write_bam_link_script(df):
    #received_not_sent_df = receivedDF[receivedDF['CRIS_Order#'].isin(received_from_earlier)][['CRIS_Order#', 'Phenotips_ID', 'Phenotips_Family_ID', 'Batch_Sent']]
    rootdir = config['rootdir']
    f = config['linkscript_fname']

    script = open(f, "w")
    script.write("#!/bin/sh\nset -e\n\n")

    for i in df.index.values.tolist():
        pbatchdir = df.loc[i]['Batch_Received']
        pbatchdir = pbatchdir.replace('BATCH0', 'BATCH')
        pid = df.loc[i]['Phenotips_ID']

        bamlink = "ln -s {}/{}/BAM/{}.bam {}/{}/BAM/\n".format(rootdir, pbatchdir, pid, rootdir, batch_name)
        bailink = "ln -s {}/{}/BAM/{}.bam.bai {}/{}/BAM/\n".format(rootdir, pbatchdir, pid, rootdir, batch_name)
        script.write(bamlink)
        script.write(bailink)
        script.write("\n")
    
    script.close()
    return()
    
def write_cumulative_ped(bsi):
    
    ped_fields = ['Phenotips_Family_ID', 
                   'Phenotips_ID', 
                   'Father_Phenotips_ID', 
                   'Mother_Phenotips_ID', 
                   'Gender', 
                   'Affected']
    
    pedDF = bsi[ped_fields]
    pedDF.Mother_Phenotips_ID.replace('^00$', '0', regex=True, inplace=True)

    # fill in missing values for Mother and Father
    pedDF.Mother_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    pedDF.Mother_Phenotips_ID.replace("", '0', inplace=True)
    pedDF.Father_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    pedDF.Father_Phenotips_ID.replace("", '0', inplace=True)
    
    pedDF.Gender.replace('M', '1', inplace=True)
    pedDF.Gender.replace('F', '2', inplace=True)

    pedDF = pedDF.drop_duplicates()
    pedDF.to_csv('cumulative_ped.txt', sep='\t', header=False, index=False)
    

####################################
# 
# Main 
#
####################################

def main():
    # Global variables
    global batch
    global batch_name
    global batch_label
    global testmode
    global config
    global pedonly

    #
    # Usage statement
    #
    parseStr = 'Reads a list of available files, performs quality control on the data,\n\
    and outputs the files needed for GRIS.\n\n\
    Usage:\n\
        csi_to_gris.py -b batch \n\n\
    Example:\n\
        csi_to_gris.py -b 7\n\
        csi_to_gris.py -b 7 -p\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--batch', required=True, type=int, help='Batch number (integer)')
    parser.add_argument('-p', '--pedfile_only', required=False, action='store_true', default=False, help='Output ped file only, do not update the tracking, masterkey or other files')
    parser.add_argument('-d', '--dir_info', required=False, type=str, default="./rawdata", help='Directory containing input CIDR csv files ("rawdata")')

    args = parser.parse_args()
    batch = args.batch
    pedonly = args.pedfile_only
    dir_info = args.dir_info

    #####################################
    #
    # Set up the variables and the log file
    #
    #####################################

    # Set up the log file
    thedate = str(datetime.datetime.now()).split()[0]
    thedate = re.sub("-","",thedate)
    global log 
    log = open('csi_to_gris' + '.log', 'a')
    log.write('\n' + str(datetime.datetime.now()) + '\n')
    log.write(' '.join(sys.argv) + '\n')
    log.write('csi_to_gris.py version ' + __version__ + '\n\n')
    log.flush()

    #####################################
    #
    # Load the Config File, and check it is complete
    #
    #####################################
    # Set Batch name (BSI value) and Batch label (output text value)
    if batch < 10: 
       batch_name = "BATCH0" + str(batch)
    else:
       batch_name = "BATCH" + str(batch)
    batch_label = "Batch " + str(batch)

    config, foundfilestext = create_config(dir_info)
    pause_for_input("\n" + foundfilestext + "Are these the correct files?\nPlease enter 'y' to continue processing the data, or 'q' to quit and correct the filenames.\n", 'y', 'q', log)

    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

        
    #####################################
    ##
    ## Import each of the reference files, do pedigree first to find the duplicate sample
    ##
    #####################################
    send_update("\nImporting information from Pedigree, Sample Mapping, Sample Key, and Manifest files...", log)

    #
    # Read the Pedigree information, find the duplicate sample id, and get the batch information
    #
#     send_update("Importing pedigree information from file: {}".format(config['pedigree']), log, True)
#     peds, theDupe = import_pedigree(config['pedigree'])
    # Subject_ID=CrisOrder#, Batch Sent
    #print("Peds {}:\n{}".format(str(peds.shape[0]), peds.head())) 

    #
    # Read the Sample Mapping information
    #
    send_update("Importing sample mapping information from file: {}".format(config['mapping']), log, True)
    mapping, theControl, collab_IDs = import_samplemapping(config['mapping'], findControl=True)
    # Subject_ID=CrisOrder#, ExomeID
    #print("Mapping Data {}:\n{}".format(str(mapping.shape[0]), mapping.head())) 

    #
    # Read the Master Sample Key information
    #
    send_update("Importing master sample key information from file: {}".format(config['samplekey']), log, True)
    samplekey = import_samplekey(config['samplekey'], theControl)
    samplecount = samplekey.shape[0]
    # Subject_ID=CrisOrder#, 6dig#
    #print("Sample Key data {}:\n{}".format(str(samplecount), samplekey.head())) 

    ########################################
    #
    # Collection information from BSI, create the output files
    #
    ########################################
    send_update("\nPulling sample data from BSI...", log)
    ## Establish a BSI Connection with user's credentials
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    ##
    ## Query for all orders included in this batch
    ##
    fields = ['CRIS Order #', 'Phenotips ID', 'Phenotips Family ID', 'Batch Sent', 'Batch Received', 'CRIS Order Status']
    sentDF = bsi_query(curl_get, url_reports, session, fields, [batch_name], 'Batch Sent')
    sentDF = sentDF[~sentDF['CRIS_Order_Status'].str.contains('Canceled')]
    sentDF = sentDF[~sentDF['CRIS_Order_Status'].str.contains('Auto Complete')]
#    sentDF.to_csv('sentDF.csv')
    #print("Sent:\n{}\n".format(sentDF.head()))

    ##
    ## Query BSI for data returned with this batch, and their family members
    ##
    # Pull data for all orders returned with Batch
    fields = ['CRIS Order #', 'Phenotips ID', 'Phenotips Family ID', 'Batch Sent', 'Batch Received', 'CRIS Order Status']
    receivedDF = bsi_query(curl_get, url_reports, session, fields, samplekey.index.tolist(), 'CRIS Order #')
#    receivedDF = receivedDF[~receivedDF['CRIS_Order_Status'].str.contains('Canceled')]
#    receivedDF = receivedDF[~receivedDF['CRIS_Order_Status'].str.contains('Auto Complete')]   
    
    
    ##
    ## Query BSI for data returned in ANY batch
    ##
    fields = ['CRIS Order #', 'Phenotips ID', 'Phenotips Family ID', 'Batch Sent', 'Batch Received', 
              'Father PhenotipsId', 'Mother PhenotipsId', 'Gender', 'CRIS Order Status', 'Affected Status']       
    all_receivedDF = bsi_query(curl_get, url_reports, session, fields, ['BATCH*'], 'Batch Received', False)     
        
#    if not pedonly:
#        receivedDF.to_csv('received.csv')
#    receivedDF.to_csv('receivedDF.csv')
    #print("Received:\n{}\n".format(receivedDF.head()))

    # Pull data for all family members in receivedDF
    fields = ['CRIS Order #', 'Batch Sent', 'Batch Received', 'Phenotips Family ID', 'Phenotips ID', 'Father PhenotipsId', 'Mother PhenotipsId', 'Gender', 'Affected Status', 'CIDR Exome ID', 'MRN', 'CRIS Order Status']
    familyDF = bsi_query(curl_get, url_reports, session, fields, receivedDF['Phenotips_Family_ID'].unique(), 'Phenotips Family ID')
    familyDF.rename(columns = {'Exome ID' : 'Exome_ID'}, inplace = True)

    ######################################
    #   
    # Creating output files (SampleTracking, Masterkey, Pedigree, Batch Info
    #
    ######################################
        
    # Write cumulative Pedigree file
    write_cumulative_ped(all_receivedDF)
        
    # Create SampleTracking
    all_family_members, sequenced_family_members, unsequenced_family_members = create_tracking(sentDF, receivedDF, familyDF, config['tracking'])

    # Print MasterKey, BatchInfo, Pedigree files
    outfiles = [config['masterkey'], config['batchinfo'], config['newpedigree']]
    write_files(outfiles, all_family_members, set(receivedDF['CRIS_Order#']), sequenced_family_members, samplekey, batch)

    send_update("\nFinished writing output files.", log)
    if not pedonly:
        write_bam_link_script(sequenced_family_members)

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    send_update("\ncsi_to_gris.py successfully completed", log)
    send_update(str(datetime.datetime.now()) + '\n', log)
    log.close()

if __name__ == '__main__':
    main()

