#Daniel Ence
#January 13, 2016
#
#Main class for reading in, filtering, and assembling short reads
#based on taxonomer inputs

#arguments to take in:
#taxonomer output
#fastqfile
#prefix for outputs/intermiediate files
#

import argparse
import sys
import os
import subprocess
from Bio import SeqIO
import re
import multiprocessing as mp
import pdb

def main():
	parser = argparse.ArgumentParser(description="Assemble reads classified to each SeqID")
	parser.add_argument("taxonomer_out",type=str,help="taxonomer output file")
	parser.add_argument("fastq",type=str,help="fastq file on which the taxonomer output is based")


if __name__ == '__main__':
	main()

