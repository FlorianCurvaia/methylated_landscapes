#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 17:30:47 2022

@author: floriancurvaia
"""
import sys


TF=str(sys.argv[1])

with open(TF+"_nomet_genonets.csv", "r") as nomet:
    with open(TF+"_met_genonets.csv", "r") as met:
        with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_comparison.csv", "w") as comp:
            lines_nomet = nomet.readlines()
            to_comp=dict()
            for row in lines_nomet:
                line=row.strip().split("\t")
                if line[0].lower() == "genotype":
                    header="Genotype\tScore_nomet\tScore_met\tDifference\tContains_CpG\n"
                    comp.write(header)
                else:
                     to_comp[line[0]]=[float(line[1])]
            
            lines_met=met.readlines()
            for row in lines_met:
                line=row.strip().split("\t")
                if line[0].lower() == "genotype":
                    pass
                elif line[0].replace("M", "C") in list(to_comp.keys()):
                    to_comp[line[0].replace("M", "C")].append(float(line[1]))
            
            for seq, scores in to_comp.items():
                if scores is not None and len(scores)==2:
                    to_write="\t".join([seq, str(scores[0]), str(scores[1]), str(scores[1]-scores[0]), str("CG" in seq)+"\n"])
                    comp.write(to_write)
            