#Daniel Ence
#July 22, 2015


from Bio import SeqIO
import os
import sys
import argparse
import re
import subprocess

#Take in a fasta file of sequences that we care about
#Take in a file of taxonomer classifier outputs
#Take a fasta file of reads to filter

parser = argparse.ArgumentParser(description="Output reads that didn't classifiy to a set of refseqs")
parser.add_argument("RefSeq",type=str,help="File of Sequences that we care about")
parser.add_argument("taxonomer_out",type=str,help="Taxonomer classifier output")
parser.add_argument("fasta",type=str,help="fasta file of reads that we want to filter")
args=parser.parse_args()

#make a hash of the the reference sequences
handle = open(args.RefSeq,'r')
refseq_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
handle.close()

#Make a hash of the taxonomer output
#If a thing didn't classifiy to one of the reference sequences, count it as unclassified 
classified_hash = {}
handle = open(args.taxonomer_out,'r')
for line in handle:
	parts = line.split("\t")
	
	if(parts[0] == "U"):
		classified_hash[parts[1]] = 1
	else:
		if(not parts[3] in refseq_dict):
			classified_hash[parts[1]] = 1			
handle.close()

#go through the fasta file and filter based on the stuff in the hash
#If something wasn't classified, output that read
fasta_file = SeqIO.parse(args.fasta,"fasta")
output_fasta_filename = os.path.basename(args.fasta).replace(".fasta",".unclassified.fasta")
output_fasta = open(output_fasta_filename,'w') 
for read in fasta_file:
	if(not read.id in classified_hash):
		SeqIO.write(read,output_fasta,"fasta")
	else:
		classified_hash.pop(read.id)
output_fasta.close()



