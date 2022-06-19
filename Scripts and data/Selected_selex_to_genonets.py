#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 13:19:55 2022

@author: floriancurvaia
"""
import sys

TF=str(sys.argv[1])

with open("Selex_analysis/"+TF+"_comparison.csv", "r") as to_gen:
    with open(TF+"_nomet_genonets.csv", "w") as gen_nomet:
        with open(TF+"_ext_genonets.csv", "w") as gen_ext:
            gen_nomet.write("Genotype"+"\t"+"Score"+"\t"+"Delta"+"\t"+"Genotypeset"+"\n")
            gen_ext.write("Genotype"+"\t"+"Score"+"\t"+"Delta"+"\t"+"Genotypeset"+"\n")
            lines=to_gen.readlines()
            seqs_scores_nomet=dict()
            seqs_scores_met=dict()
            for row in lines:
                line=row.strip().split("\t")
                if line[0]=="Genotype":
                    pass
                else:
                    seqs_scores_nomet[line[0]]=float(line[1])
                    if line[0].find("CG")!=-1:
                        seqs_scores_met[line[0]]=float(line[1])
                        seqs_scores_met[line[0].replace("CG", "MG")]=float(line[2])
                    else:
                        seqs_scores_met[line[0]]=(float(line[1])+float(line[2]))/2
            for seq_nomet, score_nomet in sorted(seqs_scores_nomet.items(), key=lambda pair:pair[1], reverse=True):
                to_write_nomet="\t".join([seq_nomet, str(score_nomet), "0.25", "set1\n"])
                gen_nomet.write(to_write_nomet)
            for seq_met, score_met in sorted(seqs_scores_met.items(), key=lambda pair:pair[1], reverse=True):
                to_write_ext="\t".join([seq_met, str(score_met), "0.25", "set1\n"])
                gen_ext.write(to_write_ext)
                
                
