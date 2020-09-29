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
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query

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
        error_text = error_text + '\nPlease enter "y" to continue or "q" to quit.\n'
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
    
    # Come back here!! Make this actually work:
    
#    collab_ped = collab_ped[~collab_ped['Subject_ID'].str.startswith("002")]
#    collab_ped['Family_ID'] = "FAM_" + collab_ped['Subject_ID'].astype(str)
#    collab_ped['Father_ID'] = "0"
#    collab_ped['Mother_ID'] = "0"
#    collab_ped['Affected'] = "Y"
#    
#    have_family = collab_ped[collab_ped['Family'] != ""]
#    fams = 
    # try using groupby() to group families and loop thru them: https://stackoverflow.com/questions/27405483/how-to-loop-over-grouped-pandas-dataframe
    
    # within the same family, if there is a # in Father/Mother, go to Individual #
        # if individual # exists in the df, then add that individual's Subject_ID to this row's Mother/Father field
        # assign Family ID as "FAM_" + min(all Subject_IDs)
            
    
    ## Each set has a sequencing duplicate that should be removed
    theDupe = find_the_duplicate(peds, "ped")

    # create a series of batch information with Subject as the index
    peds = peds[['Subject_ID', 'Investigator Column 3']]
    peds = peds[peds.Subject_ID.notnull()].set_index('Subject_ID')
    peds = peds['Investigator Column 3'].str.replace("_.*$", "", regex=True)

    if theDupe != None:
        peds = peds.drop(theDupe, errors='ignore')
        
    peds.to_csv('peds.csv')

    return(peds, theDupe)

# Import sample mapping information from csv file, export collaborator masterkey
def import_samplemapping(f, dup, findControl):
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
        collaborators = collaborators[~mapping['Subject_ID'].str.contains("NA12878")]
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

    if dup != None:
        mapping = mapping.drop(dup, errors='ignore') #if it isn't there, don't worry

    return(mapping, ctrl, collaborators)

# Import sample key information from csv file
def import_samplekey(f, dupe, ctrl):
    test_file(f)

    # import and convert to series
    samplekey = pd.read_csv(f, header=0, encoding='utf-8')
    samplekey = samplekey.set_index('Subject_ID')
    samplekey = samplekey['LIMS_SampleId']

    # remove the duplicate and control
    if dupe != None:
        samplekey = samplekey.drop(dupe)
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

    return(familyDF, seq_family_members, unseq_family_members)

# Output the three additional files
def write_files(filenames, familyDF, receivedIDs, seqFamDF, unseqFamDF):
    #filenames 0: masterkey, 1: batchinfo, 2: seqrped
    receivedDF = familyDF.loc[familyDF['CRIS_Order#'].isin(receivedIDs)]

    #seqFamDF = familyDF[familyDF['Phenotips_ID'].isin(seqIDs)]
    seqDF = pd.concat([receivedDF, seqFamDF])

    if not pedonly:
        # Write masterkey: Contains sequenced family members + those received in current batch
        master_fields = ['CRIS_Order#',
                       'Phenotips_ID', 
                       'Exome_ID',
                       'Batch_Sent', 
                       'Batch_Received']
        seqDF[master_fields].to_csv(filenames[0], sep='\t', header=True, index=False)

        # Write batchinfo
        batch_fields = ['MRN',
                        'CRIS_Order#',
                        'Phenotips_Family_ID',
                        'Phenotips_ID']
        seqDF[batch_fields].to_csv(filenames[1], sep='\t', header=True, index=False)
    
    # Write SEQR Pedigree
    ## Convert mother 00 to 0 - 00 means unknown female in BSI, 0 means unknown male, 
    ## SEQR needs both to be 0 not 00.
    #print("familyDF: {}".format(familyDF.shape))
    #print("Received: {}".format(receivedDF.shape))
    #print("UnSeq: {}".format(unseqFamDF.shape))
    #print("UnSeqFamDF: {}".format(unseqFamDF.shape))
    ped_fields = ['Phenotips_Family_ID', 
                   'Phenotips_ID', 
                   'Father_Phenotips_ID', 
                   'Mother_Phenotips_ID', 
                   'Gender', 
                   'Affected']
    seqDF.to_csv('seq_df.csv')
    
    pedDF = pd.concat([seqDF[ped_fields], unseqFamDF])
    pedDF.Mother_Phenotips_ID.replace('^00$', '0', regex=True, inplace=True)

    # fill in missing values for Mother and Father
    pedDF.Mother_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    pedDF.Mother_Phenotips_ID.replace("", '0', inplace=True)
    pedDF.Father_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    pedDF.Father_Phenotips_ID.replace("", '0', inplace=True)

    pedDF = pedDF.drop_duplicates()
    pedDF.to_csv(filenames[2], sep='\t', header=False, index=False)

    return()

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
    parser.add_argument('-d', '--dir_info', required=False, type=str, default="./rawdata/Sample_Info", help='Directory containing input CIDR csv files ("rawdata/Sample_Info/")')

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
    send_update("Importing pedigree information from file: {}".format(config['pedigree']), log, True)
    peds, theDupe = import_pedigree(config['pedigree'])
    # Subject_ID=CrisOrder#, Batch Sent
    #print("Peds {}:\n{}".format(str(peds.shape[0]), peds.head())) 

    #
    # Read the Sample Mapping information
    #
    send_update("Importing sample mapping information from file: {}".format(config['mapping']), log, True)
    mapping, theControl, collab_IDs = import_samplemapping(config['mapping'], theDupe, findControl=True)
    # Subject_ID=CrisOrder#, ExomeID
    #print("Mapping Data {}:\n{}".format(str(mapping.shape[0]), mapping.head())) 

    #
    # Read the Master Sample Key information
    #
    send_update("Importing master sample key information from file: {}".format(config['samplekey']), log, True)
    samplekey = import_samplekey(config['samplekey'], theDupe, theControl)
    samplecount = samplekey.shape[0]
    # Subject_ID=CrisOrder#, 6dig#
    #print("Sample Key data {}:\n{}".format(str(samplecount), samplekey.head())) 

    ######################################
    #
    # Compare Subject_IDs in pedigree, mapping, samplekey files
    #
    ######################################
    send_update("\n\nComparing the Subject_IDs from the pedigree, sample mapping, and sample key files...", log)
        
    # Pedigree vs Sample Key
    ped_not_key, key_not_ped = compare_keys(samplekey, peds)

    # Mapping vs Sample Key
    mapping_not_key, key_not_mapping = compare_keys(samplekey, mapping)

    # if they are different, error out until they are fixed
    # don't test ped_not_order - because these are ones released from previous batches, which we expect.
    if sum([len(ped_not_key), len(key_not_ped), len(mapping_not_key), len(key_not_mapping)]) > 0:
        send_update("Missing {} key(s) in Pedigree not Sample Key {}: ".format(str(len(ped_not_key)), ped_not_key), log)
        send_update("Missing {} key(s) in Sample Key not Pedigree {}: ".format(str(len(key_not_ped)), key_not_ped), log)
        send_update("Missing {} key(s) in Sample Mapping not Sample Key {}: ".format(str(len(mapping_not_key)), mapping_not_key), log)
        send_update("Missing {} key(s) in Sample Key not Sample Mapping {}: ".format(str(len(key_not_mapping)), key_not_mapping), log)
        
        pause_for_input("\nError: Subject_IDs in pedigree, sample mapping, sample key, and order files are inconsistent.\n" + \
                "Please enter 'y' to continue processing the data, or 'q' to quit and correct the data.\n", 'y', 'q', log)
    else:
        send_update("\nGreat News!! Subject_IDs in pedigree, sample mapping, and sample key files are consistent.\n", log)

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
    send_update("\nWriting output files...", log)
    # Create SampleTracking
    all_family_members, sequenced_family_members, unsequenced_family_members = create_tracking(sentDF, receivedDF, familyDF, config['tracking'])

    # Print MasterKey, BatchInfo, Pedigree files
    outfiles = [config['masterkey'], config['batchinfo'], config['newpedigree']]
    write_files(outfiles, all_family_members, set(receivedDF['CRIS_Order#']), sequenced_family_members, unsequenced_family_members)

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

