#Daniel Ence

import re
import os
import argparse

parser = argparse.ArgumentParser(description="Go through a gene-association file (GAF) and get a normalized count of GO terms")
parser.add_argument("gaf",type=str,help="")
args = parser.parse_args()

gaf_hash = {}	
GO_term_hash = {}	
gafFileHandle = open(args.gaf,'r')
for line in gafFileHandle:
	line_parts = re.split("\t",line)
	if(line_parts[1] in gaf_hash):
		tmp_list = gaf_hash[line_parts[1]]
		tmp_list.append(line) 
		gaf_hash[line_parts[1]] = tmp_list
	else:
		tmp_list = []
		tmp_list.append(line) 
		gaf_hash.setdefault(line_parts[1],tmp_list)
	
	if(line_parts[4] in GO_term_hash):
		GO_term_hash[line_parts[4]].append(line_parts[1])
	else:
		tmp_list = [line_parts[1]]
		GO_term_hash.setdefault(line_parts[4],tmp_list)

#get the weight for each gene_id
gene_ID_weights = {}
for gene_ID in gaf_hash:
	gene_ID_weights[gene_ID] = float(1) / float(len(gaf_hash[gene_ID]))
	#sum the weights for each GO term
for GO_term in GO_term_hash:
	curr_count = float(0)
	curr_GO_list = GO_term_hash[GO_term]
	for gene_ID in curr_GO_list:
		curr_count = float(curr_count) + float(gene_ID_weights[gene_ID])
	print GO_term + "\t" + str(curr_count)	

