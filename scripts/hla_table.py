#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:17:24 2019



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
from ncbr_bsi import read_conf, send_curl, return_bsi_info, get_bsi_name
from ncbr_huse import err_out

 # Finds batch number from the current directory 
dir_renamed = os.getcwd()
batch_number = re.sub("^.*BATCH","",dir_renamed)
batch_name = "BATCH" + str(batch_number).zfill(2)

logfile = open(dir_renamed + '/hla_log.txt', 'w')

# Get BSI connection session ID
def get_bsi_session(url, user, pw):
    curl_string = "curl -s -X POST --header 'Content-Type: application/x-www-form-urlencoded' --header 'Accept: text/plain' -d 'user_name=" + user + "&password=" + pw + "' 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'"
    sessionID = send_curl(curl_string)
#    print("Session ID: {}".format(sessionID))
#    print(sessionID.decode("utf-8"))
    
    if sessionID.decode("utf-8").find("Logon failed: The username, password, and database combination is incorrect") != -1:
            err_out("\n*** Error: login information is incorrect. ***\nQuitting.")

    #curl_string = "curl -s -X POST --header 'Content-Type: application/x-www-form-urlencoded' --header 'Accept: text/plain' -d 'user_name=" + user "' 'https://rest.bsisystems.com/api/rest/EBMS/common/logon' --digest"
    #sessionID = send_curl(curl_string)
    return(sessionID)

# Construct a query and send to BSI, return a dataframe
def bsi_query(curl, url, session, fields, batch_num, search_field, isequal=True, islike=False):
    
    fields = [get_bsi_name(f) for f in fields]

    study = "&criteria=subject.study_id%3DNIAID%20Centralized%20Sequencing"
    ## order status added as per Xi Cheng email 2/20/2019
    order_status = "&criteria=sample.field_314%3D%22Specimen%20Collected%22"
    limit = '10'
    
    # Construct the curl command
    curl += session.decode(encoding='UTF-8') + "'"
    curl += " '" + url + "?display_fields=" + "&display_fields=".join(fields) + study 
    curl += "&criteria=" + get_bsi_name(search_field) #this is the BSI name for Batch Received found in ncbr_bsi.py
    #print(curl)

    # add the "!" for not or "=@" for like
    if not isequal:
        curl += "!"
    if not islike:
        curl += "%3D" + "%3B" + batch_num
    else:
        curl += "%3D%40" + "%3B" + batch_num

    curl = curl + "&type=1'"
    #print(curl)

    # Get the data using the curl command
    data = send_curl(curl)
    # print("DATA:\n{}".format(data))

    data = data.decode('utf-8')
    data = json.loads(data)

    #print(data)

    if 'message' in data:
        if re.search("Error running report:", data['message']):
            err_out("\n*** BSI query failed to return valid results ***\nQuitting.")

    if len(data) == 0:
        err_out("BSI query failed to return valid results")
        
    # Convert the data into a dataframe
    df = pd.DataFrame(data['rows'], columns=data['headers'])
    # print("Curl results size: {}".format(df.shape))

    # Rename the columns:
    colnameDict = {'CRIS Order #': 'CRIS_Order#',
                   'Batch Sent': 'Batch_Sent',
                   'Batch Received': 'Batch_Received',
                   'Batch Ready': 'Batch_Ready',
                   'Phenotips Family ID': 'Phenotips_Family_ID',
                   'PhenotipsId': 'Phenotips_ID',
                   'Mother PhenotipsId': 'Mother_Phenotips_ID',
                   'Father PhenotipsId': 'Father_Phenotips_ID',
                   'Family Complete Status':  'Family_Complete_Status',
                   'Affected Status': 'Affected',
                   'CIDR Exome ID':'CIDR_Exome_ID',
                   'CRIS Order Status':'CRIS_Order_Status',
                   'Patient Name':'Patient_Name',
                   'Date Drawn' : 'Date_Drawn',
                   'Date Received' : 'Date_Received',
                   'GRIS Owner' : 'GRIS_Owner',
                   'Tissue Origin' : 'Tissue',
                   'Subject ID': 'MRN', 
                   'Order Date' : 'Order_Date'}
    df.rename(columns=colnameDict, inplace=True)
    df = df.replace(r'^\s*$', 'Missing', regex=True)

    return(df)

# Pull ID Dictionary from BSI
def query_BSI_data(batch_num, bsi_fields):

    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    # Returns all CIDR Exome IDs and PhenotipIDs from the specified batch
    df = bsi_query(curl_get, url_reports, session, bsi_fields, batch_num, 'Batch Received')

    return(df)

# Build HLA Table, replace Phenotips ID with CIDR Exome ID
def build_table(ID_dict, path, id_type):
    
    logfile.write('Building ' + id_type + ' HLA Table...\n')
    hla_tab = pd.DataFrame()
    cols_set = False
    colnames = []
    num_rows = ID_dict.shape[0]
 
#    Loop thru all files
    for i in range(num_rows):
        if ID_dict['is_valid_ID'][i]:
            pheno_id = ID_dict['Phenotips_ID'][i]
            fname = path + pheno_id + '/hla/R1_bestguess_G.txt'
            try:    #ID may be missing so test it
                df = pd.read_csv(fname, delimiter = "\t")
            except: #Writes out to error log and skips to next iteration
                logfile.write('Error: Phenotips ID ' + pheno_id + ' is missing from the directory.\n')
                continue
            else:   #if no error occurs, append the row 
                # concatenate ID and alleles to one row and append to HLA Table
                alleles = df['Allele'].values
                row = pd.DataFrame(np.concatenate(([ID_dict[id_type][i]], alleles), axis = 0)).T
                hla_tab = hla_tab.append(row, ignore_index = True)
                
                if not cols_set:    #if no column names, grab from first input file
                    chroms = df['Chromosome']
                    loci = df['Locus']
                    allele_cols = ["HLA-{} {}".format(locus, chrom) for chrom, locus in zip(chroms, loci)]
#                    print([id_type] + colnames)
                    colnames = [id_type] + allele_cols
                    cols_set = True
                
        
    hla_tab.reset_index(drop = True, inplace = True)    #resets row indeces to 0,1,2,3...

    coldict = dict(zip(hla_tab.columns.values.tolist(), colnames))
    hla_tab.columns = colnames
    logfile.write('\n\n')
    return hla_tab
    

def main():
    
# Usage statement and options
    parseStr = 'Creates a CSV file with Patient IDs and corresponding HLA alleles\n\
    from every patient in the entire batch. \nBy default, produces two tables: one with CIDR IDs and one with Phenotips IDs.\n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group()   #Only allowed to pick -c, -p, blank for default option
    group.add_argument('-c', '--cidr', help = 'Output one HLA table with CIDR Exome IDs', action = 'store_true')
    group.add_argument('-p', '--phenotips', help = 'Output one HLA table with Phenotips IDs', action = 'store_true')
    
    args = parser.parse_args()
  
     # Generates ID Dictionary to switch between Phenotips ID to CIDR Exome ID
    fname = 'masterkey_batch' + batch_number + '.txt'
    ID_dict = pd.read_csv(fname, sep = '\t')
    
    # Valid ID if (Batch Recieved is blank OR matches BATCH#) AND CIDR Exome ID isn't blank
    ID_dict = ID_dict.fillna(' ')
    batch_is_empty = ID_dict['Batch_Received'] == ' '
    batch_is_equal = ID_dict['Batch_Received'] == batch_name
    cidr_not_empty = ID_dict['CIDR_Exome_ID'] != ' '
    ID_dict['is_valid_ID'] = (batch_is_empty | batch_is_equal) & cidr_not_empty    #column contains all usable valid IDs

#    ID_dict = query_BSI_data(batch_name, bsi_fields) #Used to query BSI for id dictionary
    
    if args.cidr:
        hla_CIDR = build_table(ID_dict, 'HLA/', 'CIDR_Exome_ID')
        hla_CIDR.to_csv('hla_tab_cidr_batch' + batch_number + '.csv', index = False, header = True)
        print('Successfully wrote CIDR ID HLA table to ' + 'hla_tab_cidr_batch' + batch_number + '.csv')
    elif args.phenotips:
        hla_Pheno = build_table(ID_dict, 'HLA/', 'Phenotips_ID')
        hla_Pheno.to_csv('hla_tab_phenotips_batch' + batch_number + '.csv', index = False, header = True)
        print('Successfully wrote Phenotips ID HLA table to ' + 'hla_tab_phenotips_batch' + batch_number + '.csv')
    else:
        hla_CIDR = build_table(ID_dict, 'HLA/', 'CIDR_Exome_ID')
        hla_Pheno = build_table(ID_dict, 'HLA/', 'Phenotips_ID')
        hla_CIDR.to_csv('hla_tab_cidr_batch' + batch_number + '.csv', index = False, header = True)
        hla_Pheno.to_csv('hla_tab_phenotips_batch' + batch_number + '.csv', index = False, header = True)
        print('Successfully wrote HLA tables to ' + 'hla_tab_cidr_batch' + batch_number + '.csv and ' + 'hla_tab_phenotips_batch' + batch_number + '.csv')

        
    logfile.close()
   
if __name__ == '__main__':
    main()
     