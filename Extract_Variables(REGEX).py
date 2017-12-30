#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 01:17:03 2017

@author: jaspreet
"""
import os
import sys
import pandas as pd
import re

'''
ExtractVariables_RegEx(bed_file, fasta_file, mismatch) : Function that uses python Regular Expression to search a primer
in a given sequence.
Input:
    1. bed-file: Location of bed_file in format 'chr','start','end','gene','fwd_primer','rev_primer'
    2. fasta-file: Location of fasta file in format 
        >1
        ATGCATGAFGATGHAJ
        >2
        ATTTGGGTTGGTGGFFFSHDGH
    3. mismatch: Mismatch Tolerance
OUTPUT:
A tab-delimited file in the format:
        'Sequence','VariableString_Mimatch-0','VariableString_Mismatch-1','VariableString_Mismatch-2'
This uses REGEX and will provide all variables untill it reached the mismatch tolerance.
The difference in between this funtion and ExtractVariables is that ExtracVariable breaks as soon as 
it encounters a match and number of mismatch bases are <=2 but this function would output a variable for every match 
alongwith tolerance. 

re.complie(primer).finditer(sequence)- will provide (start,end) for FWD and REV primer.
If primer does not match variable==sequence[0:len(sequence)];
If fwd_primer_end="BLANK" (FWD does not match), variable = sequence[0:StartPosition_REV_Primer]
If rev_primer_end="BLANK" (REV does not match), variable = sequence[EndPosition_FWD_Primer:len(sequence)]
'''
def ExtractVariable_RegEx(bed_file, fasta_file, mismatch):
    cpd_bed = pd.read_table(bed_file, sep="\t",header=None, index_col=False, dtype=None)
    cpd_bed.columns = ['chr','start','end','gene','fwd_primer','rev_primer']
    primers = cpd_bed[['fwd_primer','rev_primer']].values.tolist()
    cpd_fasta = open(fasta_file, "r").readlines()
    variable_dict = {}
    c=0
    for line in cpd_fasta:
        line=line.rstrip()
        if re.search(">",line):
            continue
        else:
            sequence = line
            print(sequence)
            mismat=int(mismatch)
            variable_required=[]
        ## FORWARD PRIMER and REVERSE PRIMER
            prim1 = primers[c][0] 
            prim2 = primers[c][1]
            fwd = []
            rev = []
            for mm in range(mismat+1):
                if mm!=0:
                    primer1= prim1[:-(mm)]
                    res1_mm = re.compile(primer1)
                    for i_1 in res1_mm.finditer(sequence):
                        st,end = i_1.span()
                        fwd.append(end)
                    
                    primer2=prim2[:-(mm)]
                    res2_mm = re.compile(primer2)
                    for j_1 in res2_mm.finditer(sequence):
                        st,end=j_1.span()
                        rev.append(st)
        
                    if len(fwd)==0:
                        fwd.append(0)
    
                    if len(rev)==0:
                        rev.append(len(sequence))
                else:
                    res1_mm=re.compile(prim1)
                    for i_0 in res1_mm.finditer(sequence):
                        st,end=i_0.span()
                        fwd.append(end)
                        res2_mm=re.compile(prim2)
                    for j_0 in res2_mm.finditer(sequence):
                        st,end = j_0.span()
                        rev.append(st)
        
                    if len(fwd)==0:
                        fwd.append(0)
    
                    if len(rev)==0:
                        rev.append(len(sequence))
            if len(fwd)==len(rev):
                for i in range(len(fwd)):
                    variable_required.append(sequence[fwd[i]:rev[i]])
            else:
                sortfwd = sorted(fwd)
                sortrev = sorted(rev)
                if len(sortfwd)>len(sortrev):
                    fwd_1 = sortfwd[:len(sortrev)]
                fwd = fwd_1[::-1]
                for i in range(len(fwd)):
                    variable_required.append(sequence[fwd[i]:rev[i]])
        
        variable_dict[sequence]=variable_required
        c+=1
  
    df1 = pd.DataFrame.from_dict(variable_dict, orient='index')
    df1.reset_index(inplace=True)
    df1.columns = ['Sequence','VariableString_Mimatch-0','VariableString_Mismatch-1','VariableString_Mismatch-2']
    df1.to_csv("Variables_Extract(REGEX).txt", sep="\t", index=False)
        
if __name__=="__main__":
    bed_file = sys.argv[1]
    fasta_file = sys.argv[2]
    mismatch = sys.argv[3]
    ExtractVariable_RegEx(bed_file, fasta_file, mismatch)