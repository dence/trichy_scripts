#!/home/dence/applications/anaconda/bin/python
#Daniel Ence
#06/27/2015

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import test_annotate
import sys
import os
import argparse
import re
import subprocess

parser = argparse.ArgumentParser(description="annotate contigs with reference sequences")
parser.add_argument("RefSeq",type=str, help="fasta file of the reference sequences")
parser.add_argument("output",type=str, help="directory to be made and will contain the maker output directories") 
parser.add_argument("directories",nargs='+',type=str, help="list of directories containing velvet contigs")
#take a RefSeq and a bunch of assembly directories
args=parser.parse_args()

print "making directory for output " + args.output
if not os.path.exists(args.output):
	os.mkdir(args.output)

#print args.directories


#figure out which files go with which refseqs
velvet_file_dict = annotate.annotate_velvet_contigs(args.RefSeq, args.output, args.directories, 20)

annotate.batch_align(args.RefSeq, args.output, velvet_file_dict, 20)
#make a maker directory for curr refseq  
#get the velvet_contig for curr refseq from each assembly directory
#put all the contigs into a fasta file in the maker directory
#make control files with the refseq file as EST evidence and est2genome turned on 
#run maker with the contigs as the genome and the refeseq files as the EST evidence
#get the sequences with their names

