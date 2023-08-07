#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from numpy import sqrt, square
from math import isnan
import os

parser = argparse.ArgumentParser(prog="Combine profile",usage="Combine profile of circSHAPE-MaP")
parser.add_argument("-i", "--input", type=str, help="input file")
parser.add_argument("--RNA", type=str, help="Origin fa file of target RNA")
parser.add_argument("--maxbg", type=float,default=0.05,help="Max mutation ratio of background(DMSO sample),default:0.05")
parser.add_argument("--mindepth", type=int, default=1000,help="Minimal effective read depth of each samples,default:1000")
parser.add_argument("-o","--outfile", type=str, default="circRNA", help="output file name,default file is circRNA_pre-profile.txt")
args = parser.parse_args()

pd.set_option('display.max_rows',1000)
pd.set_option('display.max_columns',1000)
pd.set_option('display.width',1000)

f=pd.read_csv(args.input,sep='\t',header=0)

def Read_seq(RNA):
    with open(RNA,'r') as f:
        next(f)
        seq=next(f).strip()
    return seq
    
## RNA sequence
Origin_seq=Read_seq(args.RNA).replace("T","U")
Seq=Origin_seq+Origin_seq
seq_target=''.join(f["Sequence"].tolist())
## length of 
L=len(Origin_seq)
L_target=len(seq_target)

Pos1=Seq.index(seq_target) #target sequence start site
Pos2=Pos1+len(seq_target)   
#print Pos1,Pos2
if Pos1==-1:
    print("sequence has wrong!!!") 
    os._exit()


## Right order table
tmp1=pd.DataFrame(np.full([Pos1,f.shape[1]],np.nan), columns=f.columns.values.tolist())
New_f1=pd.concat([tmp1,f[0:L-Pos1]],axis=0).set_index(np.arange(L))

tmp2=pd.DataFrame(np.full([len(Seq)-Pos2,f.shape[1]],np.nan), columns=f.columns.values.tolist())
New_f2=pd.concat([f[L-Pos1:],tmp2],axis=0).set_index(np.arange(L))

#print New_f1
#print New_f2
print("Finish reorder table")

'''
Check sequence
'''


New_file=pd.DataFrame(np.full([L,f.shape[1]],np.nan), columns=f.columns.values.tolist())
New_file["Nucleotide"]=range(1,L+1)
New_file["Sequence"]=list(Origin_seq)
    
New_file["Modified_mutations"]=New_f1["Modified_mutations"].add(New_f2["Modified_mutations"],fill_value=0)
New_file["Modified_read_depth"]=New_f1["Modified_read_depth"].add(New_f2["Modified_read_depth"],fill_value=0)
New_file["Modified_effective_depth"]=New_f1["Modified_effective_depth"].add(New_f2["Modified_effective_depth"],fill_value=0)
#New_file["Modified_rate"]=New_file["Modified_mutations"]/New_file["Modified_effective_depth"]
        
New_file["Untreated_mutations"]=New_f1["Untreated_mutations"].add(New_f2["Untreated_mutations"],fill_value=0)
New_file["Untreated_read_depth"]=New_f1["Untreated_read_depth"].add(New_f2["Untreated_read_depth"],fill_value=0)
New_file["Untreated_effective_depth"]=New_f1["Untreated_effective_depth"].add(New_f2["Untreated_effective_depth"],fill_value=0)
#New_file["Untreated_rate"]=New_file["Untreated_mutations"]/New_file["Untreated_effective_depth"]
        
New_file["Denatured_mutations"]=New_f1["Denatured_mutations"].add(New_f2["Denatured_mutations"],fill_value=0)
New_file["Denatured_read_depth"]=New_f1["Denatured_read_depth"].add(New_f2["Denatured_read_depth"],fill_value=0)
New_file["Denatured_effective_depth"]=New_f1["Denatured_effective_depth"].add(New_f2["Denatured_effective_depth"],fill_value=0)
#New_file["Denatured_rate"]=New_file["Denatured_mutations"]/New_file["Denatured_effective_depth"]
        

'''
#### Calculate profile    
#### adjust read depth for regions where multinuc mutations prevent resolution of modifications within their length
'''
## single sample (don't use controls)
if New_file["Untreated_mutations"].isnull().any() is False and New_file["Denatured_mutations"].isnull().any() is False:
    New_file["Modified_rate"]=New_file["Modified_mutations"]/New_file["Modified_effective_depth"] # Note: requires python3 (otherwise need np.true_divide())
    New_file["Reactivity_profile"]=np.copy(New_file["Modified_rate"])
    New_file["Std_err"]=sqrt(New_file["Modified_rate"])/sqrt(New_file["Modified_effective_depth"])

# two samples (use untreated control)
elif New_file["Denatured_mutations"].isnull().any() is False:
    New_file["Modified_rate"]=New_file["Modified_mutations"]/New_file["Modified_effective_depth"]   
    New_file["Untreated_rate"]=New_file["Untreated_mutations"]/New_file["Untreated_effective_depth"]
    New_file["Reactivity_profile"]=New_file["Modified_rate"]-New_file["Untreated_rate"]
    # sqrt of sum of the squares of the individual stderrs
    stderr_M = sqrt(New_file["Modified_rate"])/sqrt(New_file["Modified_effective_depth"])
    stderr_BG= sqrt(New_file["Untreated_rate"])/sqrt(New_file["Untreated_effective_depth"])

    New_file["Std_err"]= sqrt( squart(stderr_M) + squart(stderr_BG) )
    
    
# three samples (use both controls)
else:  
    New_file["Modified_rate"]=New_file["Modified_mutations"]/New_file["Modified_effective_depth"]   
    New_file["Untreated_rate"]=New_file["Untreated_mutations"]/New_file["Untreated_effective_depth"]
    New_file["Denatured_rate"]=New_file["Denatured_mutations"]/New_file["Denatured_effective_depth"]
    New_file["Reactivity_profile"]=(New_file["Modified_rate"]-New_file["Untreated_rate"])/New_file["Denatured_rate"]
    
    stderr_M = sqrt(New_file["Modified_rate"])/sqrt(New_file["Modified_effective_depth"])
    stderr_BG= sqrt(New_file["Untreated_rate"])/sqrt(New_file["Untreated_effective_depth"])
    stderr_DC= sqrt(New_file["Denatured_rate"])/sqrt(New_file["Denatured_effective_depth"])
    
    New_file["Std_err"]=sqrt( square(stderr_M/New_file["Denatured_rate"])+\
                             square(stderr_BG/New_file["Denatured_rate"])+\
                             square( stderr_DC*(New_file["Modified_rate"]-New_file["Untreated_rate"])/square(New_file["Denatured_rate"])) )
    
#print New_file
     
# exclude nucs with high background mutation rates
# exclude nucs with low read coverage
# exclude nucs with N in the target seq
#
# depths should be a dict of numpy arrays of int
# where dict keys are: modified, untreated, denatured
New_file=New_file.replace([np.inf, -np.inf], np.nan)
New_file["HQ_profile"]=np.copy(New_file["Reactivity_profile"])
New_file["HQ_stderr"]=np.copy(New_file["Std_err"])
max_background=args.maxbg
min_depth=args.mindepth
print(min_depth)
for i in range(New_file.shape[0]):
    good_nuc = True
    if New_file["Sequence"][i].upper() not in ['A','U','G','C']:
        good_nuc = False
    if New_file["Sequence"][i].islower():
        good_nuc = False
    if New_file["Untreated_rate"].isnull().any() is False:
        if New_file["Untreated_rate"][i] > max_background:
            good_nuc = False
    if New_file["Modified_effective_depth"][i]< min_depth or New_file["Untreated_effective_depth"][i]< min_depth or New_file["Denatured_effective_depth"][i]< min_depth:
        good_nuc = False
    
    ## filter 
    if  not isnan(New_file["Reactivity_profile"][i]) and good_nuc:
        continue
    else:
        New_file.loc[i,"HQ_profile"]= np.nan
        New_file.loc[i,"HQ_stderr"] = np.nan
#print New_file


New_file["Modified_mutations"]=New_file["Modified_mutations"].astype(int)
New_file["Modified_read_depth"]=New_file["Modified_read_depth"].astype(int)
New_file["Modified_effective_depth"]=New_file["Modified_effective_depth"].astype(int)

New_file["Untreated_mutations"]=New_file["Untreated_mutations"].astype(int)
New_file["Untreated_read_depth"]=New_file["Untreated_read_depth"].astype(int)
New_file["Untreated_effective_depth"]=New_file["Untreated_effective_depth"].astype(int)

New_file["Denatured_mutations"]=New_file["Denatured_mutations"].astype(int)
New_file["Denatured_read_depth"]=New_file["Denatured_read_depth"].astype(int)
New_file["Denatured_effective_depth"]=New_file["Denatured_effective_depth"].astype(int)

New_file=New_file.replace(np.nan,"nan")
#New_file.drop(["Norm_profile"],axis=1,inplace=True)
#New_file.drop(["Norm_stderr"],axis=1,inplace=True)
New_file.to_csv(args.outfile+"_pre-profile.txt",index=False,sep='\t')  
#New_file.to_csv("test.txt",index=False,sep='\t')  
print(args.input+" write into "+args.outfile+"_pre-profile.txt")
