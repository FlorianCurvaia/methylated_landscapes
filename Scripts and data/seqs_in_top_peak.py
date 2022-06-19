#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 18:11:18 2022

@author: floriancurvaia
"""
import ast
import sys

TF=str(sys.argv[1])

deltas=[]
peaks_count_nomet=[]
l=0
with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_nomet_deltas_results.txt", "r") as res:
    lines=res.readlines()
    previous=""
    for row in lines:
        line=row.strip().split(" -")
        if l%20==0:
            deltas.append(float(line[0]))
        elif previous=="Finding peaks":
            peaks=ast.literal_eval(line[0])
            peaks_count_nomet.append(peaks["Peak 1"])
        previous=line[0]
        l+=1

peaks_count_ext=[]
with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_ext_deltas_results.txt", "r") as res:
    lines=res.readlines()
    previous=""
    for row in lines:
        line=row.strip().split(" -")
        if previous=="Finding peaks":
            peaks=ast.literal_eval(line[0])
            peaks_count_ext.append(peaks["Peak 1"])
        previous=line[0]

with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_peaks_count.txt", "w") as f:
    for d, seqs_nomet, seqs_ext in zip(deltas, peaks_count_nomet, peaks_count_ext):
        f.write("Delta: "+str(d)+"\n")
        f.write("Number of sequences nomet: "+str(len(seqs_nomet))+"\n")
        f.write("Number of sequences ext: "+str(len(seqs_ext))+"\n")
        common=len(seqs_nomet.intersection(seqs_ext))
        f.write("Number of common sequences: "+str(common)+"\n")
        unmeth=[]
        for s in seqs_ext:
            unmeth.append(s.replace("MG", "CG"))
        equiv=0
        for s in seqs_nomet:
            equiv+=unmeth.count(s)
        f.write("Number of methylated equivalent: "+str(abs(common-equiv))+"\n")