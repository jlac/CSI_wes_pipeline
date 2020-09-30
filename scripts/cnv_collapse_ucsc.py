#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:36:57 2019
Vasu Kuram
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

cnv_collapse.py
    Collapses multiple CNVs and respective fields from the same gene into one row per gene

AnnotSV Manual:
https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf

@author: Vasu Kuram
"""

import re
import sys
import datetime
import pandas as pd
import numpy as np
import argparse
import os
from ncbr_huse import test_file, err_out
from argparse import RawTextHelpFormatter

# returns all unique values if there are multiple. If not, returns the first instance of the value
def get_unique_values(repeat_values, delims, out_delim):

    vals = []

    # Split up the repeated values into individual values from the arrays of values for each gene repeat
    for value in repeat_values:
        if pd.isnull(value) == False and value != '-':
            split_vals = re.split(delims, str(value))
            vals = vals + split_vals
            
    # No values so fill in with placeholder
    if not vals:
        return None 
    
    # Only keep the unique values, separated by specified output delimiter
    output = out_delim.join(pd.Series(vals).unique())
    
    # Cleans up the output if input results in multiple spaces returned or weird commas
    if output[0] == out_delim:
        output = output[2:]
    if output[len(output) - 2:] == out_delim + ' ':
        output = output[:len(output) - 2]
            
    output = re.sub(' +', ' ', str(output))
    
    return output

# Returns the max number from values excluding -1 and '-'
def get_max(values):
    
    values = values.tolist()
    to_max = []
    
    # filter out unwanted chars
    for val in values:
        if isinstance(val, str) and ',' in val:
#            print('Its a string!')
            nums = val.split(',')
            for num in nums:
                to_max.append(num)
        elif val != '-':
            if float(val) > -1.0:
#                print(type(val))
#                print(val)
                to_max.append(val)
    if len(to_max) == 0:
        return None
    
    
    return max(to_max)

# Returns the min number from values excluding -1 and '-'
def get_min(values):

    values = values.tolist()
    to_min = []
    
    # filter out unwanted chars
    for val in values:
        if isinstance(val, str) and ',' in val:
#            print('Its a string!')
            nums = val.split(',')
            for num in nums:
                to_min.append(num)
        elif val != '-':
            if float(val) > -1.0:
#                print(type(val))
#                print(val)
                to_min.append(val)
    if len(to_min) == 0:
        return None
    
    return min(to_min)

#Deals with all GD fields: list of numbers from gd_field that correspond to each GD_ID
def get_corresp_number(cnv_tab, gene_ind, next_gene_ind, gd_field):

    #creates the dict to store GD ID as key and number as value
    keys = str(cnv_tab['GD_ID'][gene_ind]).split(';')
    vals = str(cnv_tab[gd_field][gene_ind]).split(';')
    gd_dict = dict(zip(keys, vals))
    
    #loads in other IDs and keys as long as they don't repeat
    for i in range(gene_ind + 1, next_gene_ind):
        keys = str(cnv_tab['GD_ID'][i]).split(';')
        vals = str(cnv_tab[gd_field][i]).split(';')
        key_ind = 0
        for key in keys:
            if key not in gd_dict:
                gd_dict[key] : vals[key_ind]
            key_ind = key_ind + 1
            
    return gd_dict
     
def get_keys(gd_dict):
    key_list = []
    
    for key in gd_dict.keys():
        key_list.append(str(key))
        
    return ', '.join(key_list)

def get_vals(gd_dict):
    val_list = []
    
    for val in gd_dict.values():
        val_list.append(str(val))
        
    return ', '.join(val_list)

def main():

    # Usage statement and input
    parseStr = 'Takes in CNV text file and collapses multiple rows for CNVs in the same gene into one row per gene\n\
    with a list of correpsonding CNVs. Related fields are concatenated, max/min, etc. into on row. \n\
    Output is saved into input directory to PhenotipdsID_cnv.txt\n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', required=True, nargs='?', type=str, default=None, 
                        help='Input text file containing AnnotSV IDs and multiple rows of CNVs per gene to be collapsed.')   
    parser.add_argument('-o', '--outdir', required=True, action='store', type=str, default=None, 
                        help='Directory to store the output file in.')
    
    
    args = parser.parse_args()
    infile = args.infile
    outdir = args.outdir

    # Check if file exists first
    if infile:
        test_file(infile)

    phen_ID = re.search(r'P\d{7}', infile).group()

#    outdir = infile[:re.search(r'P\d{7}.candidate', infile).span()[0]]
    logfile = open(outdir + '/cnv_' + phen_ID + '_log.txt', 'w')
    logfile.write('************ cnv_collapse.py error log: ' + phen_ID + ' ************ \n')
    cytobands = pd.read_csv('/hpcdata/dir/SCRIPTS/cytoBand.txt', delimiter = '\t', names = ['chrom', 'start', 'end', 'cytoband', 'geistain'])
    
    cnv_tab = pd.read_csv(infile, delimiter = "\t")
    hom_num = '#hom(' + phen_ID + ')'
    htz_num = '#htz(' + phen_ID + ')'
    cnv_tab.rename(columns = {hom_num:'#hom', htz_num:'#htz'}, inplace = True)
    
    fields = ['Phenotips ID', 'MergeCount', 'CytoBand'] + cnv_tab.columns.tolist()
#    fields.remove('SV length')

    final_tab = pd.DataFrame(columns = fields)  #table to be output at the end
  
    print('Analysis started:   ' + phen_ID) 
    
    # Get indeces of the first occurrence of every gene name
    all_genes = cnv_tab['Gene name'].tolist()   #list of all genes
    unique_genes = cnv_tab['Gene name'].unique()    #unique list of genes
        
    gene_indeces = []   #stores indeces (from cnv table) of the first occurrence of each gene
    for gene in unique_genes:
        gene_indeces.append(all_genes.index(gene))  # gets index of gene name in the cnv table
    
    gene_indeces.append(len(all_genes))  # allows next loop to run properly and calculate # of gene repeats
   
    for i in range(len(gene_indeces) - 1):
        gene_index = gene_indeces[i]
        next_gene_index = gene_indeces[i + 1]
        
        row = [''] * len(fields)
    
#        Make sure the CNV type is the same across all gene repeats. If not, skip iteration and exclude gene from final table 
        if cnv_tab['SV type'][gene_index:next_gene_index].nunique() == 1:
            row[fields.index('SV type')] = cnv_tab['SV type'][gene_index]
        else:
            logfile.write('Error: DEL and DUP in gene ' + all_genes[gene_index] + '... not including in output file.\n')
            continue
    
        row[fields.index('SV chrom')] = cnv_tab['SV chrom'][gene_index] #set SV chrom the same for all repeats of gene
        row[fields.index('SV start')] = cnv_tab['SV start'][gene_index] #set SV start and end as last and first repeat gene values
        row[fields.index('SV end')] = cnv_tab['SV end'][next_gene_index - 1]
        row[fields.index('SV length')] = cnv_tab['SV length'][gene_index:next_gene_index].sum()
        row[fields.index('CopyNumber')] = cnv_tab['CopyNumber'][gene_index]
        
#row[fields.index('')] = cnv_tab[''][gene_index]
        
        row[fields.index('Phenotips ID')] = phen_ID
        row[fields.index('MergeCount')] = next_gene_index - gene_index

        #query cytoband df to find where chromosome matches and start and end position match up
        chromosome = 'chr' + str(row[fields.index('SV chrom')])
        start = row[fields.index('SV start')]
        end = row[fields.index('SV end')]
        cyto_query = cytobands.query('chrom == @chromosome & start <= @start').cytoband.values
        if len(cyto_query) < 1:
            row[fields.index('CytoBand')] = 'unk'
            print('Warning: Variant not found in CytoBand data!')
        else:
            row[fields.index('CytoBand')] = chromosome + cyto_query[len(cyto_query)-1]
        
        

        #Concatenate new SV ID from updated SV info    
        row[fields.index('AnnotSV ID')] = str(row[fields.index('SV chrom')]) + '_' + str(row[fields.index('SV start')]) + '_' + str(row[fields.index('SV end')]) + '_' +row[fields.index('SV type')]                         
        row[fields.index('AnnotSV type')] = cnv_tab['AnnotSV type'][gene_index]
        row[fields.index('Gene name')] = all_genes[gene_index]
        row[fields.index('NM')] = cnv_tab['NM'][gene_index]
        row[fields.index('CDS length')] = cnv_tab['CDS length'][gene_index:next_gene_index].sum()
        row[fields.index('tx length')] = cnv_tab['tx length'][gene_index:next_gene_index].sum()
        
        #Store first and last locations if there are multiple. If not, take the first instance of the location
        if cnv_tab['location'][gene_index:next_gene_index].nunique() > 1:
            loc1 = cnv_tab['location'][gene_index]
            loc2 = cnv_tab['location'][next_gene_index - 1]
            concat_loc = ''
            
            #if the first location is past the last location, switch them when concatenating
            if re.search(r'\d+', loc1).group() > re.search(r'\d+', loc2).group():
                concat_loc = loc2[:loc2.find('-')] + loc1[loc1.find('-'):]
            else:
                concat_loc = loc1[:loc1.find('-')] + loc2[loc2.find('-'):]
            row[fields.index('location')] = concat_loc
        else:
            row[fields.index('location')] = cnv_tab['location'][gene_index]
                  
        row[fields.index('intersectStart')] = cnv_tab['intersectStart'][gene_index] #set intersect start and end as last and first repeat gene values
        row[fields.index('intersectEnd')] = cnv_tab['intersectEnd'][next_gene_index - 1]
        
        ###### DGV Gain ID Stuff ######
        row[fields.index('DGV_GAIN_IDs')] = get_unique_values(cnv_tab['DGV_GAIN_IDs'][gene_index:next_gene_index], ',', ', ')   
        row[fields.index('DGV_GAIN_n_samples_with_SV')] = get_min(cnv_tab['DGV_GAIN_n_samples_with_SV'][gene_index:next_gene_index])
        row[fields.index('DGV_GAIN_n_samples_tested')] = get_min(cnv_tab['DGV_GAIN_n_samples_tested'][gene_index:next_gene_index])
        row[fields.index('DGV_GAIN_Frequency')] = get_min(cnv_tab['DGV_GAIN_Frequency'][gene_index:next_gene_index])
        
        ###### DGV Loss ID Stuff ######
        row[fields.index('DGV_LOSS_IDs')] = get_unique_values(cnv_tab['DGV_LOSS_IDs'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('DGV_LOSS_n_samples_with_SV')] = get_min(cnv_tab['DGV_LOSS_n_samples_with_SV'][gene_index:next_gene_index])
        row[fields.index('DGV_LOSS_n_samples_tested')] = get_min(cnv_tab['DGV_LOSS_n_samples_tested'][gene_index:next_gene_index])
        row[fields.index('DGV_LOSS_Frequency')] = get_min(cnv_tab['DGV_LOSS_Frequency'][gene_index:next_gene_index])

        
        ###### GD Stuff ######
        gd_an_dict = get_corresp_number(cnv_tab, gene_index, next_gene_index, 'GD_AN')
        row[fields.index('GD_ID')] = get_keys(gd_an_dict)
        row[fields.index('GD_AN')] = get_vals(gd_an_dict)
        row[fields.index('GD_N_HET')] = get_vals(get_corresp_number(cnv_tab, gene_index, next_gene_index, 'GD_N_HET'))
        row[fields.index('GD_N_HOMALT')] = get_vals(get_corresp_number(cnv_tab, gene_index, next_gene_index, 'GD_N_HOMALT'))
        row[fields.index('GD_AF')] = get_min(cnv_tab['GD_AF'][gene_index:next_gene_index])
        row[fields.index('GD_POPMAX_AF')] = get_min(cnv_tab['GD_POPMAX_AF'][gene_index:next_gene_index])
        row[fields.index('GD_ID_others')] = get_unique_values(cnv_tab['GD_ID_others'][gene_index:next_gene_index], ';', ', ')
    

        row[fields.index('DDD_SV')] = get_unique_values(cnv_tab['DDD_SV'][gene_index:next_gene_index], ';', ', ')        
        row[fields.index('DDD_DUP_n_samples_with_SV')] = get_max(cnv_tab['DDD_DUP_n_samples_with_SV'][gene_index:next_gene_index])
        row[fields.index('DDD_DUP_Frequency')] = get_max(cnv_tab['DDD_DUP_Frequency'][gene_index:next_gene_index])
        row[fields.index('DDD_DEL_n_samples_with_SV')] = get_max(cnv_tab['DDD_DEL_n_samples_with_SV'][gene_index:next_gene_index])
        row[fields.index('DDD_DEL_Frequency')] = get_max(cnv_tab['DDD_DEL_Frequency'][gene_index:next_gene_index])
        
        row[fields.index('1000g_event')] = get_unique_values(cnv_tab['1000g_event'][gene_index:next_gene_index], ';', ', ')
        row[fields.index('1000g_AF')] = get_min(cnv_tab['1000g_AF'][gene_index:next_gene_index])
        row[fields.index('1000g_max_AF')] = get_min(cnv_tab['1000g_max_AF'][gene_index:next_gene_index])
        
        row[fields.index('IMH_ID')] = get_unique_values(cnv_tab['IMH_ID'][gene_index:next_gene_index], ';', ', ')
        row[fields.index('IMH_AF')] = get_min(cnv_tab['IMH_AF'][gene_index:next_gene_index])
        row[fields.index('IMH_ID_others')] = get_unique_values(cnv_tab['IMH_ID_others'][gene_index:next_gene_index], ';', ', ')
        row[fields.index('Mappability')] = cnv_tab['Mappability'][gene_index]
        row[fields.index('NA12878_mode_CN')] = cnv_tab['NA12878_mode_CN'][gene_index]
        row[fields.index('Repeatability_Score')] = get_min(cnv_tab['Repeatability_Score'][gene_index:next_gene_index])
        
#        These fields are mostly blank
        row[fields.index('promoters')] = get_unique_values(cnv_tab['promoters'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('dbVar_event')] = get_unique_values(cnv_tab['dbVar_event'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('dbVar_variant')] = get_unique_values(cnv_tab['dbVar_variant'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('dbVar_status')] = get_unique_values(cnv_tab['dbVar_status'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('TADcoordinates')] = get_unique_values(cnv_tab['TADcoordinates'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('ENCODEexperiments')] = get_unique_values(cnv_tab['ENCODEexperiments'][gene_index:next_gene_index], ',', ', ')

        row[fields.index('#hom')] = get_unique_values(cnv_tab['#hom'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('#htz')] = get_unique_values(cnv_tab['#htz'][gene_index:next_gene_index], ',', ', ')

        row[fields.index('GCcontent_left')] = cnv_tab['GCcontent_left'][gene_index] #set start and end as last and first repeat gene values
        row[fields.index('GCcontent_right')] = cnv_tab['GCcontent_right'][next_gene_index - 1]
        row[fields.index('Repeats_coord_left')] = cnv_tab['Repeats_coord_left'][gene_index] #set start and end as last and first repeat gene values
        row[fields.index('Repeats_coord_right')] = cnv_tab['Repeats_coord_right'][next_gene_index - 1]
        row[fields.index('Repeats_type_left')] = cnv_tab['Repeats_type_left'][gene_index] #set start and end as last and first repeat gene values
        row[fields.index('Repeats_type_right')] = cnv_tab['Repeats_type_right'][next_gene_index - 1]

        row[fields.index('ACMG')] = get_unique_values(cnv_tab['ACMG'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('HI_CGscore')] = get_min(cnv_tab['HI_CGscore'][gene_index:next_gene_index])
        row[fields.index('TriS_CGscore')] = get_unique_values(cnv_tab['TriS_CGscore'][gene_index:next_gene_index], ',', ', ')

        row[fields.index('DDD_status')] = get_unique_values(cnv_tab['DDD_status'][gene_index:next_gene_index], '/', '/')
        row[fields.index('DDD_mode')] = get_unique_values(cnv_tab['DDD_mode'][gene_index:next_gene_index], '/', '/')
        row[fields.index('DDD_consequence')] = get_unique_values(cnv_tab['DDD_consequence'][gene_index:next_gene_index], '/', ', ')
        row[fields.index('DDD_disease')] = get_unique_values(cnv_tab['DDD_disease'][gene_index:next_gene_index], '/', ', ')
        ddd_pmids = get_unique_values(cnv_tab['DDD_pmids'][gene_index:next_gene_index], '//|;|/', ', ')
        row[fields.index('DDD_pmids')] = ddd_pmids
        row[fields.index('HI_DDDpercent')] = get_min(cnv_tab['HI_DDDpercent'][gene_index:next_gene_index])

        row[fields.index('synZ_ExAC')] = get_max(cnv_tab['synZ_ExAC'][gene_index:next_gene_index])
        row[fields.index('misZ_ExAC')] = get_max(cnv_tab['misZ_ExAC'][gene_index:next_gene_index])
        row[fields.index('pLI_ExAC')] = get_max(cnv_tab['pLI_ExAC'][gene_index:next_gene_index])
        row[fields.index('delZ_ExAC')] = get_max(cnv_tab['delZ_ExAC'][gene_index:next_gene_index])
        row[fields.index('dupZ_ExAC')] = get_max(cnv_tab['dupZ_ExAC'][gene_index:next_gene_index])
        row[fields.index('cnvZ_ExAC')] = get_max(cnv_tab['cnvZ_ExAC'][gene_index:next_gene_index])
        
        row[fields.index('morbidGenes')] = get_unique_values(cnv_tab['morbidGenes'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('morbidGenesCandidates')] = get_unique_values(cnv_tab['morbidGenesCandidates'][gene_index:next_gene_index], ',', ', ')
        row[fields.index('Mim Number')] = get_unique_values(cnv_tab['Mim Number'][gene_index:next_gene_index], ';', ', ')
        row[fields.index('Phenotypes')] = get_unique_values(cnv_tab['Phenotypes'][gene_index:next_gene_index], '/|;', '; ')
        inher = get_unique_values(cnv_tab['Inheritance'][gene_index:next_gene_index], ',|//|/|;', ', ')
        
       
        # Cleans up the output if input results in multiple spaces returned
        if isinstance(inher, str):
            if inher[0] == ',':
                inher = inher[2:]
            if inher[len(inher) - 2:] == ', ':
                inher = inher[:len(inher) - 2]
                
        inher = re.sub(' +', ' ', str(inher))

        row[fields.index('Inheritance')] = inher
        row[fields.index('AnnotSV ranking')] = get_max(cnv_tab['AnnotSV ranking'][gene_index:next_gene_index])

        final_tab = final_tab.append(pd.DataFrame([row], columns = fields))
        
        
    #   Remove blank placeholders and write out the pandas df
    final_tab = final_tab.replace('-', '')
    
    outfname = outdir + '/' + phen_ID + '_cnv.txt'
    final_tab.to_csv(outfname, index = False, sep = '\t')
        
    print("{} successfully written by cnv_collapse.py".format(outfname))
    print()
    
    logfile.close()
    
if __name__ == '__main__':
    main()