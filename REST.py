#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 20:10:08 2017

@author: jaspreet
"""

import os
import sys
import urllib
import pandas as pd
from pandas import DataFrame as df
import numpy as np


'''
REST(Cancer,Genes) : Function query cBioPortal for a given TCGA-cancer data set and any number of genes simultaneously
Input: Cancer and Genes
Output: Generates "Summary_AminoAcid_Change.txt", tab-delimited file with Number of Unique Amino Acid changes wrt Genes
Detail: It query cBioPortal to :
      1. retreive Cancer Study ID for Cancer_Of_Interest(I/P) using command getCancerStudies
      2. retreive Genetic Profile ID for Cancer using command getGeneticProfiles which requires Cancer Study ID 
      as an input 
      3. retreive Mutation Data specifying AminoAcid change in corresponding Gene_Of_Interest using
      command getMutationData whihc requires Genetic Profile ID of cancer and Genes(I/P)
      4. Generates a summary text file in your working directory
'''
def REST(Cancer,Genes):        
    webservice = "http://www.cbioportal.org/webservice.do"
    url_parse = urllib.parse.urlparse(webservice) ## provides the scheme ('http'), location("www.cbioportal.org", path="/webservice.do", params, query,fragment)

    ## For retrieving all cancer study id

    query  = urllib.parse.urlencode({"cmd":"getCancerStudies"}) 
    req_urlCID = urllib.parse.urlunparse((url_parse[0],url_parse[1],url_parse[2],'',query,''))
    final_result = [i.split("\\t") for i in str(urllib.request.urlopen(req_urlCID).read()).split("\\n")]
    final_result[0][0] = final_result[0][0].split("'")[1]

    df_results = df.from_records(final_result[1:], columns=final_result[0])
    df_results.dropna(inplace=True)
    df_provsn_study = df_results[df_results['name'].str.contains("Provisional")] # TCGA DATASET
    #df_provsn_study.to_csv("Cancer_Study_ID_All_Cancers.txt", sep="\t", index=False)
    ''' Slecting id for specific Cancer'''
    Cancer_Study_ID = df_provsn_study[df_provsn_study['name'].str.contains(Cancer)]['cancer_study_id'].values.tolist()[0]

    ## For getting the Genetic Profile
    '''
    cmd=getGeneticProfiles (required)
    cancer_study_id=[cancer study ID] (required)
    '''
    query_GeneProfile = urllib.parse.urlencode({"cmd":"getGeneticProfiles","cancer_study_id":Cancer_Study_ID})
    req_urlGP = urllib.parse.urlunparse((url_parse[0],url_parse[1],url_parse[2],'',query_GeneProfile,''))
    data = [i.split("\\t") for i in str(urllib.request.urlopen(req_urlGP).read()).split("\\n")]
    data[0][0] = data[0][0].split("'")[1]
    df_resGP = df.from_records(data[1:], columns=data[0])
    df_resGP.dropna(inplace=True)
    genetic_profile_ID = df_resGP[df_resGP['genetic_profile_name'].str.contains("Mutations")]['genetic_profile_id'].values.tolist()[0]


    ## For retrieving MutationData
    ''' 
    cmd=getMutationData (required)
    genetic_profile_id= [one or more mutation profile IDs] (required)
    case_set_id= [case set ID] (optional)
    gene_list= [one or more genes, specified as HUGO Gene Symbols or Entrez Gene IDs] (required)
    '''
    #genes  = "TP53,IDH1,PIK3CA,EGFR,PTEN" ## extract genes from bed file
    query_Mutation= urllib.parse.urlencode({"cmd":"getMutationData","genetic_profile_id":genetic_profile_ID,"gene_list":Genes})
    req_MutationURL = urllib.parse.urlunparse((url_parse[0],url_parse[1],url_parse[2],'',query_Mutation,''))
    Mutation_data = [i.split("\\t") for i in str(urllib.request.urlopen(req_MutationURL).read()).split("\\n")]
    
    ## Making use of DataFrames
    df_MutationData = df.from_records(Mutation_data[2:], columns=Mutation_data[1])
    df_MutationData['Key'] = df_MutationData['chr'].astype(str)+":"+df_MutationData['start_position'].astype(str)+"-"+df_MutationData['reference_allele'].astype(str)+"-"+df_MutationData['variant_allele'].astype(str)+"|"+df_MutationData['amino_acid_change'].astype(str)
    uniqueAA_Change = []
    if len(Genes.split(","))!=1:
        for i in Genes.split(","):
            df_gene = df_MutationData[df_MutationData['gene_symbol']=="%s"%i]
            vals = df.from_records(np.unique(df_gene['Key'].values).tolist())
            lis = [i,len(vals.index)]
            uniqueAA_Change.append(lis)
    else:
        df_gene = df_MutationData[df_MutationData['gene_symbol']==Genes]
        vals = df.from_records(np.unique(df_gene['Key'].values).tolist())
        lis = [Genes,len(vals.index)]
        uniqueAA_Change.append(lis)
    df_uniq = df.from_records(uniqueAA_Change, columns=['Gene','#Unique_AminoAcid_Change'])
    df_uniq.to_csv("Summary_AminoAcid_Change.txt", sep="\t", index=False)
   
if __name__=="__main__":
    Cancer=sys.argv[1]
    Genes = sys.argv[2]
    REST(Cancer,Genes)