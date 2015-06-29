#!/home/dence/applications/anaconda/bin/python
#Daniel Ence
#06/27/2015

from Bio import SeqIO
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description="annotate contigs with reference sequences")
parser.add_argument("RefSeq",type=str, help="fasta file of the reference sequences")
parser.add_argument("output",type=str, help="directory to be made and will contain the maker output directories") 
parser.add_argument("fastas",type=str, help="txt file with a list fasta filenames. The fasta files are velvet assembled contigs")
#take a RefSeq and a bunch of assembly directories
args=parser.parse_args()

print "making directory for output " + args.output
if not os.path.exists(args.output):
	os.mkdir(args.output)

#figure out which files go with which refseqs
refseq_fasta = SeqIO.parse(args.RefSEq,"fasta") 
filesPerSeqID = {}
for refseq in refseq_fasta:
			
fastas_files = open(args.fastas,'r')
for line in fasta_files:
	line = line.rstrip("\n")

#foreach refseq
#make a maker directory for curr refseq  
#get the velvet_contig for curr refseq from each assembly directory
#put all the contigs into a fasta file in the maker directory
#make control files with the refseq file as EST evidence and est2genome turned on 
#run maker with the contigs as the genome and the refeseq files as the EST evidence
#get the sequences with their names


