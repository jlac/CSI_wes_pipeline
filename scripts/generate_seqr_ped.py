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


def generate_genrptlinks_ped(batch):
    ############ Set up the variables and the log file ############
    # Set up the log file

    thedate = str(datetime.datetime.now()).split()[0]
    thedate = re.sub("-", "", thedate)
    global log
    log = open('csi_to_gris' + '.log', 'a')
    log.write('\n' + str(datetime.datetime.now()) + '\n')
    log.write(' '.join(sys.argv) + '\n')
    log.write('csi_to_gris_hgsc.py\n\n')
    log.flush()

    # Set Batch name (BSI value) and Batch label (output text value)
    if batch < 10:
        batch_name = "BATCH0" + str(batch)
    else:
        batch_name = "BATCH" + str(batch)
    batch_label = "Batch " + str(batch)

    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    ## Establish a BSI Connection with user's credentials
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    ############ Set up output file names ############
    fnames = dict({'tracking': 'sample_tracking_summary_batch' + str(batch) + '.txt',
                   'master_key': 'masterkey_batch' + str(batch) + '.txt',
                   'ped': 'seqr_ped_batch' + str(batch) + '.txt',
                   'gen_rpt_links': 'genrptlinks_batch' + str(batch) + '.txt',
                   'link_bams': 'link_bams_batch' + str(batch) + '.sh'})

    ############ Get all active people received in this batch ############
    send_update("\nQuerying BSI for all patients received in " + batch_label + '...', log)
    fields = ['CRIS Order #', 'Phenotips ID', 'Exome ID', 'MRN', 'Phenotips Family ID', 'Batch Sent', 'Batch Received',
              'Father PhenotipsId', 'Mother PhenotipsId', 'Gender', 'CRIS Order Status', 'Affected Status',
              'Active Status']
    received = bsi_query(curl_get, url_reports, session, fields, [batch_name], 'Batch Received')

    # Remove inactive, canceled, or auto-complete orders or orders missing Exome ID
    received = received[received['Active status'] != 'Inactive']
    received = received[~received['CRIS_Order_Status'].str.contains('cancel', case=False, na=False)]
    received = received[~received['CRIS_Order_Status'].str.contains('Auto Complete', case=False, na=False)]
    received = received[~received['Exome ID'].str.contains('Auto Complete', case=False, na=False)]
    #received.to_csv('received.csv', index=False)

    ############ Get all family members of those received ############
    fam_ids = received['Phenotips_Family_ID'].unique()
    send_update("\nQuerying BSI for family members of patients received in " + batch_label + '...', log)

    # This contains all people received in the batch + all their family members (sequenced/unsequenced)
    all_fams = bsi_query(curl_get, url_reports, session, fields, fam_ids, 'Phenotips Family ID')
    all_fams = all_fams[all_fams['Active status'] != 'Inactive']

    # ped_file = all_fams[all_fams['Active status'] == 'Active']
    # ped_file = ped_file[['Phenotips_Family_ID', 'Phenotips_ID', 'Father_Phenotips_ID', 'Mother_Phenotips_ID', 'Gender', 'Affected']]
    ped_file = all_fams[['Phenotips_Family_ID', 'Phenotips_ID', 'Father_Phenotips_ID', 'Mother_Phenotips_ID', 'Gender', 'Affected']]
    genrptlinks = all_fams[['MRN', 'CRIS_Order#', 'Phenotips_Family_ID', 'Phenotips_ID']]

    # fill in missing values for Mother and Father
    ped_file.Mother_Phenotips_ID.replace('^00$', '0', regex=True, inplace=True)
    ped_file.Mother_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    ped_file.Mother_Phenotips_ID.replace("", '0', inplace=True)
    ped_file.Father_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    ped_file.Father_Phenotips_ID.replace("", '0', inplace=True)

    ped_file.drop_duplicates(inplace=True)

    return ped_file, genrptlinks, all_fams


def write_cumulative_ped():

    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    ## Establish a BSI Connection with user's credentials
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    fields = ['CRIS Order #', 'Phenotips ID', 'Phenotips Family ID', 'Batch Sent', 'Batch Received',
              'Father PhenotipsId', 'Mother PhenotipsId', 'Gender', 'CRIS Order Status', 'Affected Status']
    all_receivedDF = bsi_query(curl_get, url_reports, session, fields, ['BATCH*'], 'Batch Received', False)

    ped_fields = ['Phenotips_Family_ID',
                  'Phenotips_ID',
                  'Father_Phenotips_ID',
                  'Mother_Phenotips_ID',
                  'Gender',
                  'Affected']

    pedDF = all_receivedDF[ped_fields]
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


def main():
    pd.set_option('mode.chained_assignment', None)

    #
    # Usage statement
    #
    parseStr = 'Generates files needed for GRIS using BSI queries.\n\n\
    Usage:\n\
        csi_to_gris_bsi.py -b batch \n\n\
    Example:\n\
        csi_to_gris_bsi.py -b 7\n\
        csi_to_gris_bsi.py -b 7 -p\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--batch', required=True, type=int, help='Batch number (integer)')
    parser.add_argument('-g', '--genrptlinks', required=False, action='store_true', default=False,
                        help='Output genrptlinks file, in addition to seqr_ped')
    parser.add_argument('-c', '--cumulative', required=False, action='store_true', default=False,
                        help='Generate cumulative pedigree file in addition to seqr_ped')

    args = parser.parse_args()
    batch = args.batch
    write_genrptlinks = args.genrptlinks
    write_cumulative = args.cumulative

    ped_file, genrptlinks, all_fams = generate_genrptlinks_ped(batch)

    ############ Write out ped file ############
    ped_file.to_csv('seqr_ped_batch' + str(batch) + '.txt', header=False, index=False, sep='\t')

    if write_genrptlinks:
        genrptlinks.to_csv('genrptlinks_batch' + str(batch) + '.txt', header=False, index=False, sep='\t')
    if write_cumulative:
        write_cumulative_ped()


if __name__ == '__main__':
    main()
