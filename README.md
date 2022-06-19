# methylated_landscapes

This repository contains my python and R code to analyse genotype-phenotype landscape which contain methylated DNA sequences.

For now, the files here are not in their final version, as they need to be properly commented. Most of them are almost fully commented, but I still need to adapt the comments and finsish to comment them. 

Also, this is the first GitHub I ever created, so this might not be the cleaner GitHub you have seen, and I might reorganise it in the near future. 

Description of the main scripts:

The script sensitivity_analysis.py performs the analysis of the landscape for the desired TF as well as for the desired type of landscape (no-met or exetended). By default, analysis will be conducted for the range of delta values used in my report, but one can also specify a single value or a list of values to perform the analysis with. See the beginning of the script for example of the structure of the arguments to pass to the script.

The script treat_sensitivity_results generates the figures related to the different landscapes shown in my report. To be used, the results of both landscapes (nomet and ext) must have been generated. 

The script seqs_in_top_peak.py generates the data shown in the tables presented in the results of my report.

The script Error_spec_seq.Rmd is using the results of comparison_nomet_met.py to generate the small statistical analysis shown in my report.

SELEX_PWM_to_scores.py get scores from positions weight matrices for (Methyl-)HT-SELEX data.

Error_HT_SELEX.Rmd needs the files resulting from Comparison_SELEX.py to perform the statistical analysis shown in the report. The second chunk of the script needs data files containing all the sequences and scores for each TF, but these were too big (approx 1.5 Go in total), so the files needed to run the script starting from the third chunk are provided.

In a general way, all scripts given are actually given, as I also provided all the data files, already in the right directory, generated by each script. This allows anyone to run a specific script, without having to run the previous scirpts needed to generate the input of the given script of interest. 

All the other scripts non-mentionned in this README produced results not shown in my report. An additional description of these scripts might come in future if needed. 

Any other/previous version of the scripts displayed here can always be added on request. 

