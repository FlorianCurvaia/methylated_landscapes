#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 17:02:36 2022

@author: floriancurvaia
"""

import os 
from statistics import mean, stdev
import collections
import ast
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pl
import scipy.stats as sp
import sys
TF=str(sys.argv[1])


deltas=[]

greedy_scores_nomet=[]
random_scores_nomet=[]
peaks_count_nomet=[]
peaks_count_nomet_v2=[]
peaks_nomet=[]
stdev_greedy_nomet=[]
stdev_random_nomet=[]
all_random_nomet=[]
all_greedy_nomet=[]
all_peaks_counts=dict()
i=0
l=0
with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_nomet_deltas_results.txt", "r") as res:
    lines=res.readlines()
    previous=""
    for row in lines:
        line=row.strip().split(" -")
        if l%20==0:
            deltas.append(float(line[0]))
        elif previous=="Greedy Walk":
            greedy_scores_nomet.append(float(line[0]))
        elif previous=="Random Walk":
            random_scores_nomet.append(float(line[0]))
        elif previous=="Finding peaks":
            peaks_count_nomet.append(len(ast.literal_eval(line[0])))
            peaks_nomet.append(ast.literal_eval(line[0]))
        elif previous=="All Random Walks":
            stdev_random_nomet.append(stdev(ast.literal_eval(line[0])))
            all_random_nomet.append(ast.literal_eval(line[0]))
        elif previous=="All Greedy Walks":
            stdev_greedy_nomet.append(stdev(ast.literal_eval(line[0])))
            all_greedy_nomet.append(ast.literal_eval(line[0]))
        elif previous == "Count_greedy":
            to_n_peaks=line[0]
            n_peaks=ast.literal_eval(to_n_peaks[to_n_peaks.find("(")+1:to_n_peaks.find(")")])
            peaks_count_nomet_v2.append(len(n_peaks))
            for p, n in n_peaks.items():
                if p+"_greedy_nomet" in all_peaks_counts.keys():
                    all_peaks_counts[p+"_greedy_nomet"].append(n)
                else:
                    all_peaks_counts[p+"_greedy_nomet"]=[n]
        elif previous=="Count_random":
            to_n_peaks=line[0]
            n_peaks=ast.literal_eval(to_n_peaks[to_n_peaks.find("(")+1:to_n_peaks.find(")")])
            for p, n in n_peaks.items():
                if p+"_random_nomet" in all_peaks_counts.keys():
                    all_peaks_counts[p+"_random_nomet"].append(n)
                else:
                    all_peaks_counts[p+"_random_nomet"]=i*[0]
                    all_peaks_counts[p+"_random_nomet"].append(n)
            i+=1
        #elif previous=="Count_greedy":
        #    greedy_count.extend(ast.literal_eval(line[0]))
           
        #elif previous=="Count_random":
        #    random_count.extend(ast.literal_eval(line[0]))
    
        previous=line[0]
        l+=1

greedy_scores_ext=[]
random_scores_ext=[]
peaks_count_ext=[]
peaks_count_ext_v2=[]
peaks_ext=[]
stdev_greedy_ext=[]
stdev_random_ext=[]
all_random_ext=[]
all_greedy_ext=[]
i=0

with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_ext_deltas_results.txt", "r") as res:
    lines=res.readlines()
    previous=""
    for row in lines:
        line=row.strip().split(" -")
        if previous=="Greedy Walk":
            greedy_scores_ext.append(float(line[0]))
        elif previous=="Random Walk":
            random_scores_ext.append(float(line[0]))
        elif previous=="Finding peaks":
            peaks_count_ext.append(len(ast.literal_eval(line[0])))
            peaks_ext.append(ast.literal_eval(line[0]))
        elif previous=="All Random Walks":
            stdev_random_ext.append(stdev(ast.literal_eval(line[0])))
            all_random_ext.append(ast.literal_eval(line[0]))
        elif previous=="All Greedy Walks":
            stdev_greedy_ext.append(stdev(ast.literal_eval(line[0])))
            all_greedy_ext.append(ast.literal_eval(line[0]))
        elif previous == "Count_greedy":
            to_n_peaks=line[0]
            n_peaks=ast.literal_eval(to_n_peaks[to_n_peaks.find("(")+1:to_n_peaks.find(")")])
            peaks_count_ext_v2.append(len(n_peaks))
            for p, n in n_peaks.items():
                if p+"_greedy_ext" in all_peaks_counts.keys():
                    all_peaks_counts[p+"_greedy_ext"].append(n)
                else:
                    all_peaks_counts[p+"_greedy_ext"]=[n]
        elif previous=="Count_random":
            to_n_peaks=line[0]
            n_peaks=ast.literal_eval(to_n_peaks[to_n_peaks.find("(")+1:to_n_peaks.find(")")])
            for p, n in n_peaks.items():
                if p+"_random_ext" in all_peaks_counts.keys():
                    all_peaks_counts[p+"_random_ext"].append(n)
                else:
                    all_peaks_counts[p+"_random_ext"]=i*[0]
                    all_peaks_counts[p+"_random_ext"].append(n)
            i+=1
            
        #elif previous=="Count_greedy":
        #    greedy_count.extend(ast.literal_eval(line[0]))
           
        #elif previous=="Count_random":
        #    random_count.extend(ast.literal_eval(line[0]))
    
        previous=line[0]
"""
for sc_greedy_nomet, sc_greedy_ext, d in zip(all_greedy_nomet, all_greedy_ext, deltas):
    print(d)
    try:
        print(sp.mannwhitneyu(sc_greedy_nomet, sc_greedy_ext))
    except ValueError:
        print("not possible to compute")
for sc_random_nomet, sc_random_ext, d in zip(all_random_nomet, all_random_ext, deltas):
    
    print(d)
    try:
        print(sp.mannwhitneyu(sc_random_nomet, sc_random_ext))
    except ValueError:
        print("not possible to compute")
"""
#print(greedy_scores_ext)
#reached_peaks_greedy=collections.Counter(greedy_count)
#print(reached_peaks_greedy)
#print(random_scores_ext)
#reached_peaks_random=collections.Counter(random_count)
#print(reached_peaks_random)
with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_all_walks.csv", "w") as f:
    f.write("\t".join(["Delta","N_peaks_nomet","Greedy_nomet","Random_nomet","N_peaks_ext","Greedy_ext","Random_ext\n"]))
    for d, p_nomet, g_nomet, r_nomet, p_ext, g_ext, r_ext in zip(deltas, peaks_count_nomet_v2, all_greedy_nomet, all_random_nomet, peaks_count_ext_v2, all_greedy_ext, all_random_ext):
        for g_n, r_n, g_e, r_e in zip(g_nomet, r_nomet, g_ext, r_ext):
            f.write("\t".join([str(d),str(p_nomet),str(g_n),str(r_n),str(p_ext),str(g_e),str(r_e)+"\n"]))
        
        
plt.figure(1)
plt.clf()
#plt.plot(deltas, greedy_scores_nomet, "r", label="Greedy walk nomet")
#plt.plot(deltas, random_scores_nomet, color="blue", label="Random walk nomet")
#plt.plot(deltas, greedy_scores_ext, "g", label="Greedy walk ext")
#plt.plot(deltas, random_scores_ext, "orange", label="Random walk ext")
colors=pl.cm.viridis(np.linspace(0,1,4))
plt.axvline(0.25, ls="--", color="black", label="0.25")
plt.errorbar(deltas, greedy_scores_nomet, yerr=stdev_greedy_nomet, label="Greedy walk nomet", capsize=5, color=colors[0])
plt.errorbar(deltas, random_scores_nomet, yerr=stdev_random_nomet, label="Random walk nomet", capsize=5, color=colors[1])
plt.errorbar(deltas, greedy_scores_ext, yerr=stdev_greedy_ext, label="Greedy walk ext", capsize=5, color=colors[2])
plt.errorbar(deltas, random_scores_ext, yerr=stdev_random_ext, label="Random walk ext", capsize=5, color=colors[3])
plt.xlabel("Delta")
plt.ylabel("Sub-optimisation score")
plt.legend()
plt.savefig("Sensitivity_analysis_results/"+TF+"/"+"Plots/"+TF+"_scores.png", format="png", bbox_inches='tight', dpi=300)



plt.figure(2)
plt.clf()
plt.plot(deltas, peaks_count_nomet_v2, label="No-met")
plt.plot(deltas, peaks_count_ext_v2, label="Ext")
plt.axvline(0.25, ls="--", color="black", label="0.25")
plt.xlabel("Delta")
plt.ylabel("Number of peaks")
plt.legend()
plt.savefig("Sensitivity_analysis_results/"+TF+"/"+"Plots/"+TF+"_peaks.png", format="png", bbox_inches='tight', dpi=300)


plt.figure(3)
plt.clf()
colors_2=pl.cm.viridis(np.linspace(0,1,len(all_peaks_counts)))
plt.axvline(0.25, ls="--", color="black", label="0.25")
j=0
all_peaks_counts={k: all_peaks_counts[k] for k in sorted(all_peaks_counts)}
for peak, count in all_peaks_counts.items():
    counts=count
    if len(counts)!=len(deltas):
        to_add=len(deltas)-len(counts)
        if to_add<0:
            raise ValueError("There is more values than deltas")
        counts.extend(to_add*[0])
    plt.plot(deltas, counts, label=peak, color=colors_2[j])
    j+=1
plt.xlabel("Delta")
plt.ylabel("Number of times the peak is reached")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("Sensitivity_analysis_results/"+TF+"/"+"Plots/"+TF+"_peaks_reached.png", format="png", bbox_inches='tight', dpi=300)


with open("Sensitivity_analysis_results/"+TF+"/"+"Corr_peaks_score_"+TF+".csv", "w") as cor:
    cor.write("\t".join(["Delta","N_peaks_nomet","Greedy_nomet","Random_nomet","N_peaks_ext","Greedy_ext","Random_ext\n"]))
    for d, n_nomet, g_nomet, r_nomet, n_ext, g_ext, r_ext in zip(deltas, peaks_count_nomet_v2, greedy_scores_nomet, random_scores_nomet, peaks_count_ext_v2, greedy_scores_ext, random_scores_ext):
        to_write="\t".join([str(d), str(n_nomet), str(g_nomet), str(r_nomet), str(n_ext), str(g_ext), str(r_ext)+"\n"])
        cor.write(to_write)












