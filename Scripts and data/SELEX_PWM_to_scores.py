#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:04:37 2022

@author: floriancurvaia
"""
import numpy as np
from math import log
import sys

def cartesian(arrays, out=None): ##Source : https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    #m = n / arrays[0].size
    m = int(n / arrays[0].size) 
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
        #for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m, 1:] = out[0:m, 1:]
    return out


with open("pwm_all.csv", "r") as all_pwm:
    pwm_all = dict()
    lines=all_pwm.readlines()
    print(' '.join(lines[0].strip().split(";")).split())
    to_skip=False
    for row in lines:
        line=' '.join(row.strip().split(";")).split()
        if line[0]=="Base":
            i=0
            TF='_'.join([line[1], line[2]])
            m=line[4]
            seed=line[7]
            if TF not in pwm_all.keys():
                pwm_all[TF]=dict()
            #    to_skip=True
            #else:
            #    to_skip=False
            
        elif i==0:
            line=line[1:]
            if seed not in pwm_all[TF].keys():
                pwm_all[TF][seed]=dict()
            pwm_all[TF][seed][m]=np.zeros((4, len(line)))
            pwm_all[TF][seed][m][i]=line
            i+=1
            
        elif i == 3:
            line=line[1:]
            pwm_all[TF][seed][m][i]=line
            pwm_all[TF][seed][m] = pwm_all[TF][seed][m]/sum(pwm_all[TF][seed][m])
            
        else:
            line=line[1:]
            pwm_all[TF][seed][m][i]=line
            i+=1
        

TF=str(sys.argv[1]) #Example: "PROP1_FL"
m=str(sys.argv[2]) #Example: "Methyl-HT-SELEX"
seed=str(sys.argv[3]) #Example: "TAATTNNATTA"
pwm=pwm_all[TF][seed][m]
positions=[]
for pos in range(len(pwm[0])):
    positions.append(['A', 'C', 'G', 'T'])

if m == "Methyl-HT-SELEX":
    land_type="met"
elif m == "HT-SELEX":
    land_type="nomet"
with open(TF+"_all_"+land_type+"_genonets.csv", "w") as out: #If we account to compare different seeds for the same TF, it might be useful to write the name of the seed in the filename: with open(TF+"_"+m+"_"+seed+"_scores.csv", "w") as out:
    out.write("Genotype"+"\t"+"Score"+"\t"+"Delta"+"\t"+"Genotypeset"+"\n")
    seqs=cartesian(positions)
    scores=np.prod(cartesian(np.transpose(pwm)), axis=1)
    for sequence, score in sorted(zip(seqs, scores), key=lambda pair: pair[1], reverse=True):
        sc=score
        seq=''.join(sequence)
        if m=="Methyl-HT-SELEX":
            seq=seq.replace("CG", "MG")
        if score==0:
            sc=1.e-140
        out.write(seq+"\t"+str(log(sc))+"\t"+"0"+"\t"+"set1"+"\n")
    out.close()
    
        
        
