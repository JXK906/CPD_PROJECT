#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 01:24:16 2017

@author: jaspreet
"""

import sys
import os
import pandas as pd
import re

''' 
A. GetPrimerPositions(sequence, mismatch, primer) : Function matches the primer sequence and provide a list as 
Output : [start_position_ofprimer,end_position_ofprimer, mismatch_count]

Inputs required:
    1. sequence : Fasta sequence containing variables and primers
    2. mismatch : Mismatch Tolerance
    3. primer   : Primer sequence 
'''

def GetPrimerPositions(sequence, mismatch, primer):
        global result, wp, count, shortsequence, pseudosequence
        mm=int(mismatch)
        wp = (len(sequence)-len(primer))+mm # iteration param
        result = []
        for i in range(wp):
             if (len(primer)+i)<=(len(sequence)):
                 count=0
                 shortsequence = sequence[i:len(primer)+i]
                 if shortsequence!=primer:
                     for j in range(len(primer)):
                         if primer[j]!=shortsequence[j]:
                             count+=1
                             if count>mm:
                                 break
                             else:
                                 continue
                     if count<=mm:
                        result.append([i,len(primer)+i,count])
                        break
                 else:
                    result.append([i,len(primer)+i, count])
             else:
                pseudosequence = sequence+"X"
                #print(pseudosequence)
                count=0
                shortsequence = pseudosequence[i:len(primer)+i]
                #print(shortsequence)
                if shortsequence!=primer:
                    for j in range(len(primer)):
                        if primer[j]!=shortsequence[j]:
                            count+=1
                            if count>mm:
                                #print("No Match found")
                                break
                            else:
                                continue
                    if count<=mm:
                        result.append([i,len(primer)+i, count])
                        break
                else:
                    result.append([i,len(primer)+i, count])
        return result

'''
B. ExtractVariables(bed_file,fasta_file, mismatch) : Function provides a tab-delimited file showing 
Original_Sequence_with_Primers and Extracted_Sequence/Variable_after_Removing_Primers.

If the primers do not match it considers the sequence as a variable itself. If only fwd primer matches it will extract
sequence[endpos_fwd_primer:len(sequence)] AND if only reverse primer matches it will extract sequence[0:startpos_rev_primer]
and if both matched the function extracts sequence[endpos_fwd_primer:startpos_rev_primer]

endpos_fwd_primer =  Last iterated position of FWD_primer -- #of Mismatched Bases in FWD_Primer 

This function breaks as soon as it encounters a match and number of mismatch bases are <=2

Input following:
    1. bed_file : Location of Tab-delimited file with 'chr','start','end','gene','fwd_primer','rev_primer' specified.
    2. fasta_file: Location of Fasta format
    3. mismatch:  Mismatch Tolerance as required by user
'''

# Open Bed file and extract primer sequences in a dictionary ["Fwd_Primer":"Rev_primer"]
def ExtractVariables(bed_file,fasta_file, mismatch):
    cpd_bed = pd.read_table(bed_file, sep="\t",header=None, index_col=False, dtype=None)
    cpd_bed.columns = ['chr','start','end','gene','fwd_primer','rev_primer']
    primers = cpd_bed[['fwd_primer','rev_primer']].values.tolist()
    cpd_fasta = open(fasta_file, "r").readlines()
    mm=int(mismatch)
    variable_dict = {}
    c=0
    for line in cpd_fasta:
        line=line.rstrip()
        if re.search(">",line):
            continue
        else:
            sequence = line
            variable_required=[]
            prim1 = primers[c][0]
            prim2 = primers[c][1]
            fwd_result = GetPrimerPositions(sequence, mm, prim1)
            rev_result = GetPrimerPositions(sequence, mm, prim2)
            if len(fwd_result)==1 and len(rev_result)==1:
                FWD_bases_not_matched = fwd_result[0][2] ## mismatch between primer and sequence
                FWD_last_iterated_position= fwd_result[0][1] ## len(primer)+i 
                FWD_Primer_Endposition = FWD_last_iterated_position-FWD_bases_not_matched
                REV_Primer_Startposition= rev_result[0][0] ## len(primer)+i 
                variable_required.append(sequence[FWD_Primer_Endposition:REV_Primer_Startposition])
            
            elif len(fwd_result)==0 and len(rev_result)==0:
                variable_required.append(sequence)
                
            elif len(fwd_result)>=1 and len(rev_result)==0:
                for i in fwd_result:
                    FWD_bases_not_matched = i[2] ## mismatch between primer and sequence
                    FWD_last_iterated_position= i[1] ## len(primer)+i 
                    FWD_Primer_Endposition = FWD_last_iterated_position-FWD_bases_not_matched
                    variable_required.append(sequence[FWD_Primer_Endposition:])
            elif len(fwd_result)==0 and len(rev_result)>=1:
                for i in rev_result:
                    REV_Primer_Endposition = i[0]
                    variable_required.append(sequence[REV_Primer_Endposition:])
        variable_dict[sequence]=variable_required
        c+=1
    
    outfile = open("Variables.txt",'w')
    outfile.write("#"+"\t"+"Fasta_Sequence"+"\t"+"Variables_after_removing_primers"+"\n")
    d=1
    for key, values in variable_dict.items():
        outfile.write(str(d)+"\t"+str(key)+"\t"+str(values)+"\n")
        d+=1
    outfile.close()


#ExtractVariables("/Users/jaspreet/CPD/cpd.bed","/Users/jaspreet/CPD/cpd.fasta", "2")

if __name__ == "__main__":
    bed_file = sys.argv[1]
    fasta_file = sys.argv[2]
    mismatch = sys.argv[3]
    ExtractVariables(bed_file,fasta_file,mismatch)
    
