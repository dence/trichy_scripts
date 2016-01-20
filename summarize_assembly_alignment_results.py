#Daniel Ence 
#August 7, 2015

from Bio import SeqIO

import sys
import os
import argparse
import subprocess

#take in the est2genome fasta file
#count the number of sequences, find out how long they are


parser = argparse.ArgumentParser(description="summary the number and length of est2genome hits, the input for the alignments")
parser.add_argument("RefSeq",type=str,help="fasta file of all the sequences you made alignemnts for")
parser.add_argument("prefix",type=str,help="prefix for the output files that will be created")
parser.add_argument("est2genome_files",nargs='+',type=str,help="list of est2genome fasta files")

args = parser.parse_args()


prefix = args.prefix
fasta_handle = SeqIO.parse(open(args.RefSeq,'r'),"fasta")
length_hash = {}
presence_hash = {}
sample_hash = {}
for refseq in fasta_handle:
	
	length_hash.setdefault(refseq.id,0)
	presence_hash.setdefault(refseq.id,0)

	ref_length = len(refseq.seq)	
	seq_count = 0
	est2genome_file = refseq.id + ".est2genome.fasta"
	est2genome_handle = ""
	if(not os.path.exists(est2genome_file)):
		continue

	est2genome_handle = SeqIO.parse(open(est2genome_file,'r'),"fasta")
		
	curr_refseq_lengths = {}
	curr_refseq_presence = {}
	for record in est2genome_handle:
		#There will be a copy of the the refseq sequence in 
		#the est2genome files, we don't want to get nubmers from that one
		if(record.id != refseq.id):
			sample_id = record.id.split("_")[0]
			curr_refseq_lengths[sample_id] = '{:.2}'.format( float(len(record.seq)) / float(ref_length))
			curr_refseq_presence[sample_id] = 1
			sample_hash.setdefault(sample_id, 1)
	
	length_hash[refseq.id] = curr_refseq_lengths
	presence_hash[refseq.id] = curr_refseq_presence

#now to print it out
header_line = "SeqID\t"
sample_list = sample_hash.keys()
for sample in sample_list:
	header_line = header_line + sample + "\t"									

print "presence/absence of sequences"
presence_file = prefix + ".presence.summary.txt"
presence_handle = open(presence_file,'w')
print header_line + "\n"
presence_handle.write(header_line)
for refseq in presence_hash:
	line = refseq + "\t"
	results = presence_hash[refseq]

	if(results == 0):
		results = {}
		for sample in sample_list:
			results[sample] = 0

	for sample in sample_list:
		if( not results.has_key(sample)):
			line = line + "0\t"
		else:
			line = line + str(results[sample]) + "\t"
	print line
	presence_handle.write(line + "\n")	
presence_handle.close()

print "Percent length of reference seq"
percent_length_file = prefix + ".length.summary.txt"
percent_length_handle = open(percent_length_file,'w')
print header_line
percent_length_handle.write(header_line + "\n")
for refseq in length_hash:
	line = refseq + "\t"
	results = length_hash[refseq]

	if(results == 0):
		results = {}
		for sample in sample_list:
			results[sample] = 0

	for sample in sample_list:
		if(not results.has_key(sample)):
			line = line + "0\t"
		else:
			line = line + str(results[sample]) + "\t"
	print line
	percent_length_handle.write(line + "\n")
percent_length_handle.close()
