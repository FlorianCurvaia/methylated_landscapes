#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 10:40:02 2022

@author: floriancurvaia
"""
import sys

TF=str(sys.argv[1])

with open(TF+"_best_nomet_genonets.csv", "r") as nomet:
    with open(TF+"_best_met_genonets.csv", "r") as met:
        with open("Selex_analysis/"+TF+"_comparison.csv", "w") as comp:
            
            lines_nomet = nomet.readlines()
            lines_met=met.readlines()
            sequences_met=set()
            sequences_nomet=set()
            i=0
            for row_nomet in lines_nomet:#, row_met in zip(lines_nomet, lines_met):
                line_nomet=row_nomet.strip().split("\t")
                #line_met=row_met.strip().split("\t")
                if line_nomet[0].lower() == "genotype":
                    header="Genotype\tScore_nomet\tScore_met\tDifference\tContains_CpG\n"
                    comp.write(header)
                elif i<1000:
                    seq_nomet=line_nomet[0]
                    #seq_met=line_met[0].replace("MG", "CG")
                    sequences_nomet.add(seq_nomet)
                    #sequences_met.add(seq_met)
                    i+=1
            #sequences=sequences_nomet.intersection(sequences_met)
            to_comp=dict()
            for row_nomet in lines_nomet:
                line_nomet=row_nomet.strip().split("\t")
                if line_nomet[0].lower() == "genotype":
                    pass
                elif line_nomet[0] in sequences_nomet:
                    to_comp[line_nomet[0]]=[float(line_nomet[1])]
            for row_met in lines_met:
                line_met=row_met.strip().split("\t")
                if line_met[0].lower() == "genotype":
                    pass
                else:
                    seq=line_met[0].replace("MG", "CG")
                    if seq in sequences_nomet:
                        to_comp[seq].append(float(line_met[1]))
            
            for seq, scores in to_comp.items():
                if scores is not None and len(scores)==2:
                    to_write="\t".join([seq, str(scores[0]), str(scores[1]), str(scores[1]-scores[0]), str("CG" in seq)+"\n"])
                    comp.write(to_write)