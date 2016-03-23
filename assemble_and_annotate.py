#Daniel Ence
#

import argparse
import sys
import os
#from TaxonomerRunner import TaxonomerRunner
from FastqFilter import FastqFilter
from SlowFastqFilter import SlowFastqFilter

parser = argparse.ArgumentParser(description="Classify reads and perform reference-free assembly of gene models.")
parser.add_argument("fastq",type=str,help="fastq file of reads")
#add protein nucleotide switch for classifer
parser.add_argument("output",type=str,help="Name of working directory and prefix for output files")
parser.add_argument("RefSeq",type=str,help="Fastq file with the sequences that we care about")
parser.add_argument("taxonomer",type=str,help="Output file from taxonomer")


args = parser.parse_args()

#First, make a working directory and
#Expand this out to a logger class??
log_file = ""
working_dir = ""
if not os.path.exists(args.output):
	os.mkdir(args.output)
	working_dir = args.output
log_file = "./" + working_dir + "/" + args.output + ".log"
log_handle = open(log_file,'w')
log_handle.write("Started run in this directory\t" + working_dir)

fastqFilter = SlowFastqFilter()
#fastqFilter = SlowFastqFilter()
fastqFilter.loadFastqFile(args.fastq)
fastqFilter.loadTaxonomerOutput(args.taxonomer,args.RefSeq)
fastqFilter.filterFastq(args.RefSeq,working_dir)
newFastqList = fastqFilter.getNewFastqList()


 
