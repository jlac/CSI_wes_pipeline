#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:07:30 2018
ncbr_huse.py
    Set of functions supporting the FNL NCBR work
    
"""

__author__ = 'Susan Huse'
__version__ = '1.0.0'
__copyright__ = 'none'

#import csv
import sys
import os
import re
#import datetime
import subprocess
import pandas as pd
import urllib
import json
# import requests
# from requests.auth import HTTPDigestAuth
from ncbr_huse import send_update, err_out, pause_for_input

####################################
# 
# Functions for connecting to BSI
#
####################################
def read_conf(cnf):
    if not os.path.isfile(cnf):
        # SUE ADD Log to err_out
        err_out("Error: unable to locate BSI authentication file:  " + cnf)

    with open(cnf, 'r') as f:
        x = f.readlines()

    user = x[0].rstrip()
    pw = x[1].rstrip()

    return(urllib.parse.quote(user, safe=''), urllib.parse.quote(pw, safe=''))

# Get the BSI credentials, URLs, cnf, etc
def return_bsi_info():
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "
    return(cnf, url_session, url_reports, curl_get)

# Submit the curl string to system to execute
def send_curl(curl_string):
    proc = subprocess.Popen([curl_string], stdout = subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    if err is not None:
        err_out("Errored out attempting to establish session with BSI: .{}".format(err))
    return(out)

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

# Get the table and field code names within BSI for query construction
def get_bsi_name(infield):
    fieldDict = {'CRIS Order #' : 'sample.field_274',
        'Phenotips ID' : 'sample.field_252',
        'Phenotips ID Subject' : 'subject_131.field_173',
        'Phenotips Family ID' : 'subject_131.field_170',
        'Seqr ID' : 'subject_131.field_254',
        'Batch Sent' : 'sample.field_323',
        'Batch Sent Subject': 'subject_131.field_194',
        'Batch Received' : 'sample.field_324',
        'Batch Received Subject' : 'subject_131.field_195',
        'Batch Ready' : 'sample.field_340',
        'Instructive Case' : 'subject_131.field_220', 
        'Instructive Case Comments' : 'subject_131.field_221',
        'Vendor' : 'sample.field_306',
        
        'Father PhenotipsId' : 'subject_131.field_161',
        'Mother PhenotipsId' : 'subject_131.field_167',
        'Father Phenotips ID' : 'subject_131.field_161',
        'Mother Phenotips ID' : 'subject_131.field_167',
        'Family Complete Status': 'subject_131.field_188',
        'Family MRN' : 'subject_131.field_157',
        'Father MRN' : 'subject_131.field_160',
        'Mother MRN' : 'subject_131.field_166',

        'Adopted' : 'subject_131.field_149',
        'Relationship' : 'subject_131.field_182',
        'Affected Status' : 'subject_131.field_150',
        'Affected' : 'subject_131.field_150',
        'Active Status' : 'subject_131.field_189',
        'Archive' : 'subject_131.field_216',
        'CRIS Report Genes' : 'subject_131.field_253',

        'CMA' : 'subject_131.field_203',
        'Exome ID' : 'sample.field_337',
        'CIDR Exome ID' : 'sample.field_337',
        'DLM LIS Number' : 'sample.field_336',
        'CRIS Order Status' : 'sample.field_314',
        'MRN' : 'sample.subject_id',
        'Date Drawn' : 'sample.date_drawn',
        'GRIS Owner' : 'subject_131.field_191',
        'Patient Name' : 'sample.field_322',
        'Patient Name Subject' : 'subject_131.field_169',
        'Date Received' : 'vial.date_received',
        'Gender' : 'subject_131.field_163',
        'Sex' : 'sample.sex',
        'Race' : 'sample.field_326',
        'Ethnicity' : 'subject_131.field_155',
        'Age' : 'sample.field_208',
        'Tissue Origin' : 'vial.tissue_origin',
        'Proband' : 'subject_131.field_171', 
        'Order Date' : 'sample.field_297',
        'Date of Birth' : 'subject_131.field_152',

        'First drafter' : 'subject_131.field_204',
        'Date report drafted' : 'subject_131.field_205',
        'Date team was notified' : 'subject_131.field_208',
        'Date documented in CRIMSON' : 'subject_131.field_214',
        'Report Status' : 'subject_131.field_206',
        'Discloser' : 'subject_131.field_207',
        'Date disclosed to patient' : 'subject_131.field_211',
        'Date of 1st contact' : 'subject_131.field_209',
        'Date of 2nd contact' : 'subject_131.field_210',
        'Date of CRIS upload' : 'subject_131.field_212',
        'Date data in Illumina' : 'subject_131.field_259',
        'Date of Enrollment' : 'subject_131.field_196',
        'Date data returned' : 'subject_131.field_257',
        'Date uploaded to seqr' : 'subject_131.field_258',

        'Box' : 'location.box',
        'Row' : 'vial_location.row',
        'Col' : 'vial_location.col'
        }
    return(fieldDict[infield])

# Construct a query and send to BSI, return a dataframe
def bsi_query(curl, url, session, fields, theIDs, search_field, isequal=True, islike=False):
    fields = [get_bsi_name(f) for f in fields]

    study = "&criteria=subject.study_id%3DNIAID%20Centralized%20Sequencing"
    ## order status added as per Xi Cheng email 2/20/2019
    order_status = "&criteria=sample.field_314%3D%22Specimen%20Collected%22"
    limit = '10'
    
    # replace spaces in the IDs with "%20"
    theIDs = [re.sub(" ", "%20", x) for x in theIDs]

    # Construct the curl command
    curl += session.decode(encoding='UTF-8') + "'"
    ## order status added as per Xi Cheng email 2/20/2019
    curl += " '" + url + "?display_fields=" + "&display_fields=".join(fields) + study 
    #curl += " '" + url + "?display_fields=" + "&display_fields=".join(fields) + study + order_status
    curl += "&criteria=" + get_bsi_name(search_field)
    #print(curl)

    # add the "!" for not or "=@" for like
    if not isequal:
        curl += "!"
    if not islike:
        curl += "%3D" + "%3B".join(theIDs)
    else:
        curl += "%3D%40" + "%3B".join(theIDs)

    curl = curl + "&type=1'"
#    print(curl)

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

    return(df)