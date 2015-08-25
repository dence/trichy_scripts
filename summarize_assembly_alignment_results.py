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
parser.add_argument("est2genome_files",nargs=+,type=str,helpo="list of est2genome fasta files")


args = parser.parse_args()



fasta_handle = SeqIO.parse(open(args.RefSeq,'r'),"fasta")
length_hash = []
presence_hash = []
for refseq in fasta_handle:
	
	est2genome_file = refseq.id + ".est2genome.fasta"
	est2genome_handle = SeqIO.parse(open(args.est2genome_file,'r'),"fasta")
	
	ref_length = length(refseq.Seq)	
	seq_count = 0
	for record in est2genome_handle:
		#There will be a copy of the the refseq sequence in 
		#the est2genome files, we don't want to get nubmers from that one
		if(record.id != refseq.id):
			
											
					
		


