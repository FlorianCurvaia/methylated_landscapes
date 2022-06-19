#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 28 15:06:22 2022

@author: floriancurvaia
"""

import sys
import igraph as ig #Python library to handle graphs
from itertools import permutations #Function that will help us to create an iterator object containing all the possible pairs of sequences
import numpy as np #Library to handle arrays. Provide exp function that avoids a "division by zero" error, outputs a simple warning instead, in exp(-x) when x is too close to 0
import matplotlib.pyplot as plt # Library to do various type of plots. Will be used for horizontal bar plots
from statistics import mean #We import the function to compute the mean of a list
import random as rdm #Library to generat pseudo-random numbers, as well as pseudo-random sampling
import time #Library to compute the time need for each main part of the script. Might be useful to target parts to time-optimize
import collections #Library that allows to get the count of each unique items in a list
import warnings #To be able to raise warnings during the execution of the script
from math import exp #!!Unused!! exp function. has been replaced by np.exp to avoid error
import networkx as nx #Another python library to handle graph. Useful to find the dominant network in case of disconnected graphs. (needed because function of igraph to do so is currently not working)
import os
import ast


TF=str(sys.argv[1]) #Inline argument 1: name of the transcription factor, and if needed followed by "_" and the sample number. 
#Example: CTCF_1 stands for TF CTCF sample 1
#When treating a HT-SELEX file, one can also add "_" + either "FL" or "eDBD" to differentiate for a same TF if the data were genrated using the full length TF (FL) or only the DNA binding domain (eDBD)
land_type=str(sys.argv[2]) # Inline argument 2: the type of landscape desired. Either nomet, for a landscape containing only non-methylated sequences, or ext for the extended landscape containing both methylated and non-methylated sequences
if land_type!= "nomet" and land_type != "ext":
    raise ValueError("Invalid landscape type was given. The valid types are nomet and ext. "+land_type+" was given")
# It is possible to give a third argument, which is a list of the delta one wants to test. The argument must be entered as a python list between " ". Example: "[0, 0.1, 0.2, 0.3, 0.4, 0.5]"
if len(sys.argv) < 4: #If no deltas were specified
    deltas=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
else: #If a list of deltas was given
    deltas=ast.literal_eval(sys.argv[3])


os.makedirs(os.getcwd()+"/Sensitivity_analysis_results/"+TF+"/Plots", exist_ok=True)

class mseq:
    """ 
    creation of a class to handle sequences of methylDNA alphabet as network's node, which is standard DNA alphabet + M for methyl-cytosine
    """
    def __init__(self, sequence, score, node_id, d):
        self.seq=sequence #Store the DNA/Methyl-DNA sequence of the node
        self.score=score #Store the affinity of the sequence
        self.id=node_id #Store the id of the node (by default 0 is the id of the best binder)
        self.in_neighbors=set() #Neighbors for which an edge emerge from neighbor to self
        self.out_neighbors=set() #Neighbors for which an edge emerge from self to neighbor
        self.identical_seqs=set() #Other mseq object for which their sequence is identical to self.seq
        self.is_peak=False #Is the sequence of the mseq object part of a peak in the landscape
        self.in_deviation=set() #Other mseq object for which their score is at most delta greater or smaller
        self.indels=set() #neighbours mseq thanks to indel
        self.is_valley=False # Is the sequence of self part of a valley ? !! Unused for now!!
        self.sup_neighbors=set() #Neighbors of self that have a score strictly greater to self.score
        self.delta=d
    def distance(self, potential_neighbor):
        """
    
        Parameters
        ----------
        potential_neighbor : other object of the mseq class
        
            Compute the number of mutation needed to go from self.seq to potential_neighbor.seq
            

        Returns
        -------
        Mutational distance from self to neighbor

        """
        if len(self.seq)!=len(potential_neighbor.seq): # The sequences we want to compare, as well as the sequences that we want in our landscape must be of same length 
            raise ValueError("Sequences must be of same lenght")
        h=0 #Initializing a modified Hamming distance to 0 
        for i in range(len(self.seq)): #Parsing through self.seq
            nuc_1=self.seq[i]
            nuc_2=potential_neighbor.seq[i]
            if nuc_1!=nuc_2: #If the two sequences differs at one position we need to increase h
                if nuc_2=="M": #If in sequence 2 there is a methylated C (denoted M) we need to consider different mutationnal path
                    if nuc_1=="C": #To go from C to M (methylation), we consider it being only one step
                        h+=1
                    else: #For all others mutations, we consider it as two mutations, first to C, then to M
                        h+=2
                else: #If the difference between the sequences doesn't involve a M, just add 1
                    h+=1
        if h==0: #If the distance is 0, it means that both sequences are the same
            self.identical_seqs.add(potential_neighbor.id)
            potential_neighbor.identical_seqs.add(self.id)
        return(h)
    def is_methylated(self, other_seq):
        if self.seq.replace("CG", "MG")==other_seq.seq or other_seq.seq.replace("CG", "MG")==self.seq:
            return True
        else:
            return False
    def deviation(self, other_seq):
        """
        

        Parameters
        ----------
        other_seq : mseq object
            
            find if two mseq object are in a 0.25 range within each other, with respect to their scores.

        Returns
        -------
        None.

        """
        
        if abs(self.score-other_seq.score)<=self.delta:
                self.in_deviation.add(other_seq.id)
                other_seq.in_deviation.add(self.id)
                
    def edge(self, other_seq):
        """
        

        Parameters
        ----------
        other_seq : other object of the mseq class
            
            Compute the edges between self and other_seq

        Returns
        -------
        list of all edges between self and other_seq
        None if no edge

        """
        dist=self.distance(other_seq) #Get the mutational distance from self.seq to other_swe
        if dist>1 and not self.is_indel(other_seq) and not self.is_methylated(other_seq): #If the distance is more than 1 and both sequences are not indel of each other, there is no edges between the two nodes
            return(None)
        elif dist==0: #if the sequences are identical, which should not happen, as identical sequences are considered are collapsed into a single sequence before computing edges
            self.out_neighbors.add(other_seq.id)
            other_seq.in_neighbors.add(self.id)
            self.in_neighbors.add(other_seq.id)
            other_seq.out_neighbors.add(self.id)
            return([(self.id, other_seq.id), (other_seq.id, self.id)])
        else:
            self.deviation(other_seq)
            if other_seq.id in self.in_deviation: #Test if the score of the two sequences are different, i.e. if the absolute value of their difference is greater than 0.25
                #If the absolute value of their difference is equal or smaller than 0.25, add a double arrowed edge between these nodes
                self.out_neighbors.add(other_seq.id)
                other_seq.in_neighbors.add(self.id)
                self.in_neighbors.add(other_seq.id)
                other_seq.out_neighbors.add(self.id)
                if self.score-other_seq.score<0: #Test which sequence has a score which is greater than the other
                    self.sup_neighbors.add(other_seq.id)
                elif self.score-other_seq.score>0:
                    other_seq.sup_neighbors.add(self.id)
                return([(self.id, other_seq.id), (other_seq.id, self.id)])
            diff=self.score-other_seq.score #Get the difference in the scores when it is respectively greater/smaller than 0.25/-0.25
            if diff >0:
                # Create the edge in a single direction: from other_seq to self
                self.in_neighbors.add(other_seq.id)
                other_seq.out_neighbors.add(self.id)
                other_seq.sup_neighbors.add(self.id)
                return([(other_seq.id, self.id)])
            else:
                # Create the edge in a single direction: from self to other_seq
                self.sup_neighbors.add(other_seq.id)
                self.out_neighbors.add(other_seq.id)
                other_seq.in_neighbors.add(self.id)
                return([(self.id, other_seq.id)])
            
    def is_indel(self, other_seq): #Check the sequences of two mseq object are indels of size 1
        if self.distance(other_seq)>1:
            c=0
            d=0
            for i in range(len(self.seq)-1):
                if self.seq[i]==other_seq.seq[i+1]:
                    c+=1
            for i in range(1, len(self.seq)):
                if self.seq[i]==other_seq.seq[i-1]:
                    d+=1
            if c==len(self.seq)-1 or d==len(self.seq):
                self.indels.add(other_seq.id)
                other_seq.indels.add(self.id)
                return(True)
            else:
                return(False)
        else:
            return(False)
        
            
with open("Sensitivity_analysis_results/"+TF+"/"+TF+"_"+land_type+"_deltas_results.txt", "w") as d_res:
    for delta in deltas:
        d_res.write(str(delta)+"\n")
        test=mseq("TATMGGAC", 0.34, 100000, 0.25) #Useless, just to test that the class mseq is working and show how it works
        nodes=[] #List that will contains the mseq object for each unique sequence of the input file
        scores=[] #List that will contains the scores for each unique sequence of the input file
        sequences=[] #List that will contains the sequence for each unique sequence of the input file
        is_part_of_peak=[] #List of boolean that will contains True if sequence is part of a peak(actually part of the summit of the peak), else it contains False; for each unique sequence of the input file
        #A file with the following name must be in the directory: TF_ + land_type + "_genonets.csv"
        with open(TF+"_"+land_type+"_genonets.csv", "r") as gen: #Input file in the genonets' input file format with the column order as follows: Genotype, Score, Delta, Genotypeset. 
        #The sequences' alphabet must be A, T, G, C, M. M is used only for methylated sequences to denotes methylated cytosine in a CpG dinucleotide
        #Also, the sequences order in the file should be in descending order of (relative) binding affinity (or a proxy for it)
            lines=gen.readlines()
            i=0 #Will serve as mseq.id. Will grow with the number of lines, thus the closer to 0 is the mseq.id, the better the mseq.seq binds the TF
            for row in lines:
                line=row.strip().split("\t")
                if line[0]=="Genotype":
                    pass
                else:
                    nodes.append(mseq(line[0], float(line[1]), i, delta)) #TODO: find a way take the delta from the input file, instead of modifying it directly in the script. Explanations: We define each mseq object that we will need
                    #We fill the list that we'll need with the relevant informations
                    scores.append(float(line[1]))
                    sequences.append(line[0])
                    is_part_of_peak.append(False) #At the beginning all entries are False, they'll be switched to True for sequences belonging to a peak later during computation of the peaks
                    i+=1
                    
        start_time_0 = time.time()           
        for two_nodes in permutations(nodes, r=2): #Compute de distances between the nodes of all possible pairs of nodes to detect and get identical sequences. Idea to optimize: get the neighbours at this moment, but might be problematic when needing to remove duplicate sequences 
            two_nodes[0].distance(two_nodes[1])
        
        d_res.write("Computing distances --- %s seconds ---\n" % (time.time() - start_time_0))
        
        
        start_time_1 = time.time()
        same_seqs=dict() #Dictionnary that will store all sets of identical sequences
        for node in nodes: #Parse through all the mseq object in the purpose to reduce idenctical sequences to a single sequence
            if len(node.identical_seqs)!=0: #Check if there is identical sequences
                try: #Try to see if the sequence is already in same_seqs[key], if yes, add the new sequence
                    same_seqs[node.seq]= same_seqs[node.seq].union(node.identical_seqs)
                except KeyError: #if the sequence is not already in same_seqs[key], create the entry
                    same_seqs[node.seq]=node.identical_seqs
                    
        
        
        to_remove=[] #List with the index of the nodes that we will need to remove because of duplicity
        for IDs in list(same_seqs.values()): 
            to_keep=sorted(IDs)[0] #We keep one copy of each object present multiple times
            scores[to_keep]=mean([scores[i] for i in IDs])
            nodes[to_keep].score=scores[to_keep]
            to_remove.extend(sorted(IDs)[1:])
            
        if len(to_remove)!=len(set(to_remove)): #Make sure that each mseq object to remove is present only once
            raise ValueError("some identical sequences have been considered twice or more")
        
        to_remove=sorted(set(to_remove), reverse=True) #We sort the list of objects to remove by descending order to keep to avoid to change the index of object that will be remove later
        for index in to_remove: #We delete de cuplicate objects
            del nodes[index]
            del scores[index]
            del sequences[index]
            del is_part_of_peak[index]
            
        nodes.sort(key=lambda x: x.score, reverse=True)
        for node in nodes: #We reset the nodes id to avoid out of range indices
            node.id=nodes.index(node)
        scores=sorted(scores, reverse=True)
        d_res.write("Clearing graph --- %s seconds ---\n" % (time.time() - start_time_1))
        
        start_time_2 = time.time()
        edges=set()
        
        for two_nodes in permutations(nodes, r=2): #We compute all the edges of the, now clean, graph
            ed=two_nodes[0].edge(two_nodes[1])
            if ed != None:
                for e in ed:
                    edges.add(e)
        
        # Check if there is no isolated nodes:
        is_isolated=[]
        for node_to_check in nodes:
            is_isolated.append(len(node_to_check.out_neighbors)+len(node_to_check.in_neighbors))
        if 0 in is_isolated:
            warnings.warn("Warning...........There is some isolated nodes")
        d_res.write("Finding edges --- %s seconds ---\n" % (time.time() - start_time_2))
        
        start_time_3 = time.time()
        treated_nodes=set()
        peaks=dict()
        peaks_score=dict()
        valleys=dict() #Unused, but can be used to store valleys the same way as peaks is being used
        j=1 #count the number of peaks
        k=1 #Unused, I don't remember what it is for
        for node in nodes:
            if node.id in treated_nodes: #To avoid treating the same node multiple times
                pass
            else:
                if (node.out_neighbors.issubset(treated_nodes) and len(node.out_neighbors)!=0) or len(node.out_neighbors.intersection(treated_nodes))!=0: #If a node has in its out_neighbors a node already treated, it cannot be a peak
                    treated_nodes.add(node.id)
                else:
                    to_peak=set()
                    treated_nodes.add(node.id)
                    really_a_peak=True #Will be set to False if node is not a peak, i.e. if it's part of a plateau
                    peaks_neighbors=set() #Will collect the sequences that are in the neighbourhood of the peak if they have a score higher than a sequence of the peak from which they are a neighbour. If not empty it means the peak is actually a plateau
                    height=node.score #Maximal height of the peak, because nodes are treated in descending order of their score, so the node with the highest score of a peak is always treated first
                    to_investigate=node.in_neighbors #Nodes to investigate to know if they are part of the peak or not
                    while len(to_investigate)!=0: #While there is still nodes that might belong to the peak, keep going
                        to_investigate_further=set() #This will be the nodes that when all the nodes in to_investigate will be treated that we'll need to treat next
                        for i in to_investigate:
                            if abs(height-scores[i])<=delta: #Test if the node coule be part of the peak
                                to_peak.add(nodes[i])
                                
                                if min(nodes[i].out_neighbors)<node.id:
                                    peaks_neighbors.add(i)
                                #for nei in nodes[i].out_neighbors-set([node.id])-node.in_neighbors: #Parse the out_neighbours of the nodes which is in the peak, to see if some of them might also be in the peak
                                #   if scores[nei]>scores[i]+0.25 or nodes[nei].is_peak==True: 
                                #       peaks_neighbors.add(nei)
                                to_further=[j for j in nodes[i].in_neighbors if abs(height-scores[j])<=delta and j not in treated_nodes]
                                to_investigate_further=to_investigate_further.union(set(to_further)) #Add the neighbours to investigate_further
                            treated_nodes.add(nodes[i].id) #Node has been treated
                        to_investigate=to_investigate_further #Investigate all the neighbours of the nodes that were previously investigated
                    if len(peaks_neighbors)!=0: #If peaks_neighbours is not empty, then it's not a peak but a plateau
                        really_a_peak=False
                    else:#if really_a_peak: #If the detected peak is really a peak, we fill the peaks and peaks_score dictionaries and set mseq.is_peak=True for all sequences in the peak
                        peaks["Peak "+str(j)]=set()
                        peaks["Peak "+str(j)].add(node.seq)
                        peaks_score["Peak "+str(j)]=set()
                        peaks_score["Peak "+str(j)].add(node.score)
                        node.is_peak=True
                        is_part_of_peak[node.id]=True
                        for p in to_peak:
                            peaks["Peak "+str(j)].add(p.seq)
                            peaks_score["Peak "+str(j)].add(p.score)
                            p.is_peak=True
                            is_part_of_peak[p.id]=True
                        j+=1
                    
        max_of_each_peak=dict() #Get the maximal score of each peak
        for p, seqs_of_peak in peaks_score.items():
            max_of_each_peak[max(seqs_of_peak)]=p
            
        d_res.write("Finding peaks --- %s seconds ---\n" % (time.time() - start_time_3))
        d_res.write(str(peaks)+"\n")
        d_res.write(str(peaks_score)+"\n")
        
        seq_to_max_of_peak=dict() #Look up table to look for the max value of the peak to which a sequence belong
        for seqs, sco in zip(peaks.values(), peaks_score.values()):
            max_of_peak=max(sco)
            for seq in seqs:
                seq_to_max_of_peak[seq]=max_of_peak
        
        landscape=ig.Graph(n=len(nodes), edges=list(edges), directed=True, vertex_attrs={"score":scores, "sequence":sequences, "is_peak":is_part_of_peak}) #Create graph object that contains the genotype network
        d_res.write("Is the graph connected ? : "+str(landscape.is_connected(mode="weak"))+"\n") #Check if graph is weakly connected (meaning for directed graph that the corresponding undirected graph is fully connected)
        if not landscape.is_connected(mode="weak"): #If the graph is not weakly connected we use the dominant network as a graph and redefine all the lists needed to create the graph accordingly
            list_edges=landscape.get_edgelist()
            G=nx.MultiDiGraph(list_edges)
            l_sub_g=max(nx.weakly_connected_components(G), key=len) #Get the dominant network
            l_scores=[]
            l_nodes=[]
            l_sequences=[]
            l_is_part_of_peak=[]
            for i in l_sub_g:
                l_scores.append(scores[i])
                l_nodes.append(nodes[i])
                l_sequences.append(sequences[i])
                l_is_part_of_peak.append(is_part_of_peak[i])
            scores=l_scores
            nodes=l_nodes
            sequences=l_sequences
            is_part_of_peak=l_is_part_of_peak
            for node in nodes:
                node.id=nodes.index(node)
                node.out_neighbors=set()
                node.in_neighbors=set()
                node.sup_neighbors=set()
            edges=set()
        
            for two_nodes in permutations(nodes, r=2):
                ed=two_nodes[0].edge(two_nodes[1])
                if ed != None:
                    for e in ed:
                        edges.add(e)
            landscape=ig.Graph(n=len(nodes), edges=list(edges), directed=True, vertex_attrs={"score":scores, "sequence":sequences, "is_peak":is_part_of_peak})
            
        
        ## Greedy walk:
        start_time_4 = time.time()
        rdm.seed(0) #Random seed for reproducibility
        max_affinity=scores[0] #Highest score (= score of the highest peak)
        range_max_min_peak=max_affinity-min(seq_to_max_of_peak.values()) #!!Unused!! Compute the difference between the height of the highest peak and the height of the lowest peak
        reached_peak_greedy=[] #Will contain the peak reached from each nodes in the greedy walk
        diff_peak_greedy=[] #Will contain the difference in height between heighest peak and the peak reached for eahc starting node of the greedy walk
        global_maximum_accessibility_greedy=[]
        paths_greedy=[] #All paths that have been used in the greedy walk
        
        for rep in range(100):
            visited_nodes=[] #Nodes alreday visited during greedy walk
            for node in nodes:
                if node.id in visited_nodes: #If the node has already been visited, no need to recompute a path from it as in a greedy walk there should be no major changes
                    pass
                else:
                    actual_node=node
                    actual_path=[]
                    actual_path.append(node.id)
                    visited_nodes.append(actual_node.id)
                    new_nodes=[] #All the nodes that have not been visited previously but that will be during this path
                    new_nodes.append(actual_node.id)
                    while not actual_node.is_peak: #While a pek has not been reached, keep going in the path 
                        neighbours_id=list(actual_node.out_neighbors) #Get the neighbours to which it's posssible to go from the current node 
                        neighbours_scores=[scores[i] for i in neighbours_id] #Get the scores of the reachable nodes
                        score_max=max(neighbours_scores) #Get the heighest of the scores of the reachable nodes
                        max_scores=[] #Will contain the set of all sequences that have a similar score to the heighest score
                        n=0
                        for sc in neighbours_scores: #Find all sequences that have a similar score to the heighest score
                            if abs(score_max-sc<=delta):
                                max_scores.append(neighbours_id[n])
                            n+=1
                        if len(max_scores)>1: #If there is more than one possible node to select (meaning that there is similar scores to the highest score)
                            next_node=nodes[rdm.sample(max_scores, k=1)[0]] #Randomly choose a node among the highest scores
                            if next_node.id not in visited_nodes: #If the node has not been visited previously, add it to new_nodes
                                actual_path.append(next_node.id) #Should this be out of the if condition ?
                                visited_nodes.append(next_node.id)
                                new_nodes.append(next_node.id)
                            actual_node=next_node #Make one step in the path
                        else: #If there is only one possible node for the highest score among neighbours
                            next_node=nodes[max_scores[0]]
                            if next_node.id not in visited_nodes: #If the node has not been visited previously, add it to new_nodes
                                actual_path.append(next_node.id) #Should this be out of the if condition ?
                                visited_nodes.append(next_node.id)
                                new_nodes.append(next_node.id)
                            actual_node=next_node #Make one step in the path
                    paths_greedy.append(actual_path)
                    diff_peak_greedy.extend([max_affinity-seq_to_max_of_peak[actual_node.seq]]*len(set(new_nodes))) #We add n times the difference between highest peak and reached peak with n being the number of nodes in the path that hadn't been visited before. If we want to standardize the difference in peak's height with respect to the range among with peaks' heights are distributed: /range_max_min_peak)
                    reached_peak_greedy.extend([max_of_each_peak[seq_to_max_of_peak[actual_node.seq]]]*len(set(new_nodes))) #We add n times the reached peak with n being the number of nodes in the path that hadn't been visited before.
            
            global_maximum_accessibility_greedy.append(sum(diff_peak_greedy)/len(diff_peak_greedy)) #We compute the sub-optimisation score for the greedy walk through the landscape
        d_res.write("Greedy Walk --- %s seconds ---\n" % (time.time() - start_time_4))
        d_res.write(str(sum(global_maximum_accessibility_greedy)/len(global_maximum_accessibility_greedy))+"\n")
        d_res.write("Count_greedy\n")
        #d_res.write(str(reached_peak_greedy)+"\n")
        peak_count_greedy=collections.Counter(reached_peak_greedy) #We take the number of times each peak has been reached 
        
        d_res.write(str(peak_count_greedy)+"\n")
        d_res.write("All Greedy Walks ---\n")
        d_res.write(str(global_maximum_accessibility_greedy)+"\n")
        
        #Ploting the horizontal barplots of the number of times where, and frequencies at which, each peak was reached.
        plt.figure(1)
        plt.clf()
        n_per_peak=list(peak_count_greedy.values())
        width=list(map(lambda x: x/len(nodes), n_per_peak))
        n_count=[str(round(pct, ndigits=2))+"%, n = "+str(val) for val, pct in zip(n_per_peak, width)]
        labels = list(peak_count_greedy.keys())
        y_pos=np.arange(len(labels))
        plt.barh(y_pos, width, height=0.5)
        plt.yticks(y_pos, labels)
        plt.subplots_adjust(left=0.35)
        plt.xlim(0,100)
        for i in range(len(y_pos)):
            plt.text(x=width[i]+5, y=y_pos[i], s=n_count[i])
        plt.title("Number of time each peak is reached in 100 greedy walks")
        plt.savefig("Sensitivity_analysis_results/"+TF+"/Plots/Barh_greedy_"+land_type+"_"+str(delta)+".png", format="png", bbox_inches='tight', dpi=300)
        
        #TODO: Takes a lot of time to run for CTCF extended sample 1. Might be due to the fact that it's too unlikely to move when 2 is too small. Should we standardize s with respect to all the pfix of the enighbours to gain time ?
        ## Random Walk
        def methyl_diff(seq1, seq2):##TODO: look that it doesn't return 10 for 10 #Function to weight the probability of picking a neighbour depending on if it differs by methylation or not 
            for n_1, n_2 in zip(seq1, seq2):
                #Maybe could fix the todo : #if seq1 in seq2.indels: return 1 #But for that, need to pass mseq object to the function, and not the seq
                if n_1 != n_2:
                    change=n_1+n_2
                    if change=="CM" or change=="MC":
                        return 10
                    elif change=="MT":
                        return 6
                    else:
                        return 1
            
            raise ValueError("There is still identical sequences in the graph")
        
        def p_fix(node_1, node_2): #Function to compute the probability of fixation of the allele picked among the neighbours. For a population of 100
            s=node_2.score-node_1.score
            div=1-np.exp(-4*s*100)
            if div==0:
                div=1.e-15
            return((1-np.exp(-2*s))/div)
        
        node_to_peak_path={k:[] for k in range(len(nodes))} #Store for each node the path that was done #To fix, only store the results of the last random walx because it's out of the for loop
        node_to_peak_reached={k:[] for k in range(len(nodes))} #Store for each node the peak that was reached  // ???really??? ->#To fix, only store the results of the last random walk because it's out of the for loop
        node_to_peak_diff={k:[] for k in range(len(nodes))} #Store for each node the difference between the highest peak and the peak reached #To fix, only store the results of the last random walx because it's out of the for loop
        global_maximum_accessibility_random=[] #We store the 100 sub-optimisation scores for the random walk through the landscape
        
        no_go_nodes=set() #Nodes where to not go because they are "False" peak in a plateau. To avoid oscillations
        reached_peak_random=[] #Peak that was reached
        start_time_5 = time.time()
        for rep in range(100):
            diff_peak_random=[] #Difference in height between the highest peak and the peak reached
            visited_nodes=set() #Nodes that were visited and
            #paths_random=[] #Paths that were used during the walk
            for node in nodes:
                if node.id in visited_nodes: #Avoid computig the path from the same node more than once
                    pass
                else:
                    actual_node=node
                    #actual_path=[] #Path from starting node to peak
                    #actual_path.append(node.id) #Starting node of the path
                    visited_nodes.add(actual_node.id) #We visited this node
                    if len(actual_node.sup_neighbors)==0 and actual_node.is_peak!=True: #If a node has no neighbours with a highest score than its score, but is not par of a peak, then it's a "false" peak in a plateau
                        no_go_nodes.add(actual_node.id)
                    while not actual_node.is_peak:
                        to_next=False
                        #to_pass=False
                        neighbours_id=list(actual_node.in_neighbors.union(actual_node.out_neighbors))
                        #for n in neighbours_id:
                        #    if n in no_go_nodes:
                        #        neighbours_id.remove(n)
                        
                            
                        #if len(neighbours_id)==0 and len(actual_node.out_neighbors)!=0: #To avoid being stucked in a flat part which is not a peak because of the deviation of 0.25
                        #    neighbours_id=list(actual_node.out_neighbors)
                        neighbours_pfix={i:p_fix(actual_node, nodes[i]) for i in neighbours_id}
                        methyl_weight_prob=[methyl_diff(actual_node.seq, nodes[i].seq) for i in neighbours_id]
                        #no_meth_change=[a*b for a,b in zip(neighbours_scores, methyl_weight_prob)]
                        #no_met_tot=sum(no_meth_change)
                        if actual_node.id in no_go_nodes:
                            next_node=nodes[rdm.choices(neighbours_id, weights=methyl_weight_prob, k=1)[0]]
                            to_next=True
                            
                        while to_next==False:
                            next_node=nodes[rdm.choices(neighbours_id, weights=methyl_weight_prob, k=1)[0]]
                            """
                            while to_pass==False and (len(next_node.sup_neighbors)==0 and next_node.is_peak!=True):
                                if len(neighbours_id)==0:
                                    neighbours_id=list(set(list(actual_node.in_neighbors)+list(actual_node.out_neighbors)))
                                    methyl_weight_prob=[methyl_diff(actual_node.seq, nodes[i].seq) for i in neighbours_id]
                                    next_node=nodes[rdm.choices(neighbours_id, weights=methyl_weight_prob, k=1)[0]]
                                    to_pass=True
                                    to_next=True
                                elif set(no_go_nodes).union(set(actual_node.sup_neighbors))==set(no_go_nodes):
                                    next_node=nodes[rdm.choices(neighbours_id, weights=methyl_weight_prob, k=1)[0]]
                                    to_pass=True
                                    to_next=True
                                else:
                                    pos_to_remove=neighbours_id.index(next_node.id)
                                    no_go_nodes.append(next_node.id)
                                    neighbours_id.remove(next_node.id)
                                    del methyl_weight_prob[pos_to_remove]
                                    next_node=nodes[rdm.choices(neighbours_id, weights=methyl_weight_prob, k=1)[0]]
                                    """
                            if len(next_node.sup_neighbors)==0 and next_node.is_peak!=True:
                                to_next=True
                                no_go_nodes.add(next_node.id)
                            
                            elif rdm.uniform(0, 1) <= neighbours_pfix[next_node.id]:
                                to_next=True
                            
                        #actual_path.append(next_node.id)
                        #visited_nodes.append(next_node.id)
                        actual_node=next_node
                        
                    #node_to_peak_path[node.id].append(actual_path)
                    #paths_random.append(actual_path)
                    #node_to_peak_diff[node.id].append((max_affinity-seq_to_max_of_peak[actual_node.seq])/range_max_min_peak)
                    reached_height=seq_to_max_of_peak[actual_node.seq]
                    diff_peak_random.append((max_affinity-reached_height))#/range_max_min_peak)
                    #node_to_peak_reached[node.id].append(seq_to_max_of_peak[actual_node.seq])
                    reached_peak_random.append(max_of_each_peak[reached_height])
            
            global_maximum_accessibility_random.append(sum(diff_peak_random)/len(diff_peak_random))
        d_res.write("Random Walk --- %s seconds ---\n" % (time.time() - start_time_5))
        d_res.write(str(sum(global_maximum_accessibility_random)/len(global_maximum_accessibility_random))+"\n")
        d_res.write("Count_random\n")
        #d_res.write(str(reached_peak_random)+"\n")
        peak_count_random=collections.Counter(reached_peak_random)
        d_res.write(str(peak_count_random)+"\n")
        d_res.write("All Random Walks ---\n")
        d_res.write(str(global_maximum_accessibility_random)+"\n")
        
        
        plt.figure(2)
        plt.clf()
        n_per_peak=list(peak_count_random.values())
        width=list(map(lambda x: x/len(nodes), n_per_peak))
        n_count=[str(round(pct, ndigits=2))+"%, n = "+str(val) for val, pct in zip(n_per_peak, width)]
        labels = list(peak_count_random.keys())
        y_pos=np.arange(len(labels))
        plt.barh(y_pos, width, height=0.5)
        plt.yticks(y_pos, labels)
        plt.subplots_adjust(left=0.35)
        plt.xlim(0,100)
        for i in range(len(y_pos)):
            plt.text(x=width[i]+5, y=y_pos[i], s=n_count[i])
        plt.title("Number of time each peak is reached in 100 random walks")
        plt.savefig("Sensitivity_analysis_results/"+TF+"/Plots/Barh_random_"+land_type+"_"+str(delta)+".png", format="png", bbox_inches='tight', dpi=300)
    
