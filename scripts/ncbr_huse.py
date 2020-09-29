#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:07:30 2018

ncbr_huse.py
    Set of functions supporting the FNL NCBR work
    
"""

__author__ = 'Susan Huse'
__version__ = '1.0.0'
__copyright__ = 'none'

import csv
import sys
import os
import re
import datetime
import subprocess
import pandas as pd
 
####################################
# 
# Functions 
#
####################################

#
# Run any command and send notice to log file.  if dorun=F, just notify for testing
#
def run_cmd(theCommand, fn, dorun):
    # Run the given command, x is a list of command parameters
    # if dorun = false, just send out the notification of the next python step
    print(" ".join(theCommand) + "\n")
    fn.write(" ".join(theCommand) + "\n")
    fn.flush()
    if dorun:
        subprocess.check_call(theCommand, stderr=fn)

#
# Run an OS command and notify log file
#
def run_os_cmd(theCommand, fn, dorun):
    # Run the given command, x is a list of command parameters
    # if dorun = false, just send out the notification of the next python step
    theCommandStr = " ".join(theCommand)
    print(theCommandStr + "\n")
    fn.write(theCommandStr + "\n")
    fn.flush()
    if dorun:
        os.system(theCommandStr)

#
# gunzip a file
#
def un_gzip(fname, logfn):
    # If rerunning and the previous step has already been compressed,
    # need to uncompress it before you can rerun the command
    gzip_name = fname + '.gz'
    if (not os.path.isfile(fname)) and os.path.isfile(gzip_name):
        ungz_cmd = ['gunzip', gzip_name]
        run_cmd(ungz_cmd, logfn, True)


#
# Read ~/.my.cnf and connect to an SQL database
#
def con_db(host_name, db_name, port_number):
    with open(os.path.expanduser("~/.my.cnf")) as f:
        ## [client]
        ## user="user_name"
        ## password="password"

        cnfdata = f.readlines()
        cnfdata = [x.strip() for x in cnfdata] 
        cnfdata = [x.replace('"', '') for x in cnfdata]
        user_name = cnfdata[1].replace('user=', '')
        password = cnfdata[2].replace('password=', '')
    
    db = MySQLdb.connect(host=host_name, db=db_name, port=port_number, user=user_name, passwd=password)
    return(db)

#
# Print updates to screen and log file
#
def send_update(updateStr, log=None, quiet=False):
    if not quiet:
        print(updateStr)

    if log is not None:
        log.write(updateStr + '\n')
    return(0)

#
# Log error message and exit
#
def err_out(errMsg, log = None):
    if log is not None:
        log.write(errMsg)

    sys.exit(errMsg)

#
# Test validity of input file or directory
#
def test_file(f, log=None):
    if type(f) != str:
        f = f.name
    if not os.path.isfile(f):
        err_out("Error: unable to locate input file:  " + f, None)
    else:
        return(True)

def test_dir(d, log=None):
    if not os.path.isdir(d):
        err_out("Error: unable to locate input directory:  " + d, None)
    else:
        return(True)

#
# Pause for user to be ready to continue, use contkey=None to get any input
#
def pause_for_input(txt, contkey='y', quitkey='q', log=None):

    # tally the number of tries
    ans_cnt = 0 

    # loop for the user to enter input, give them a few tries
    while True:

        # wait for the input
        ans = input(txt)

        # check if valid continue or quit keys
        if quitkey is not None and ans == quitkey:
            err_out("User elected to quit.  Exiting...\n", log)

        # if none, just return the input
        if contkey is None:
            return(ans)

        # if there is a contkey, then be sure it is correctly typed
        elif ans == contkey:
            return(ans)

        else:
        # give them additional help and increment the answer count
            reminder = "Note: only {} to continue and {} to quit are valid options.\nPlease try again.\n".format(contkey, quitkey)
            if ans_cnt == 0:
                txt = "\n" + txt + "\n" + reminder

            # Otherwise 3 strikes and exit from the loop
            if ans_cnt == 2:
                err_out("User failed to continue ({}) or quit ({}) three times in a row.  Exiting...".format(contkey, quitkey), log)

            ans_cnt = ans_cnt + 1 


#
# Count sequences in a fasta file
#
def fasta_count(fastaFile):
    seqcount = 0
    for line in open(fastaFile, "r"):
         if re.match(">", line):
               seqcount += 1
    return(seqcount)

#
# Count sequences in a fasta file
#
def fasta_list(fastaFile):
    seqs = []
    for line in open(fastaFile, "r"):
         if re.match(">", line):
               seqs.append(re.sub(">", "", line.rstrip()))
    return(seqs)

