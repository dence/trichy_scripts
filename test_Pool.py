#from multiprocessing import Pool
#from multiprocessing import Process
#import os

#def f(x):
#	return x*x
#
#if __name__ == '__main__':
#	p=Pool(5)
#	print (p.map(f,[1,2,3,4,5,5,6,6,5,5,5,7,7,10,100]))

#def info(title):
#	print title 
#	print 'module name:',__name__
#	if hasattr(os, 'getppid'): # only on Unix
#		print 'parent process:',os.getppid()
#	print 'process id:', os.getpid()

#def f(name):
#	info('function f')
#	print 'hello',name

#if __name__ == '__main__':
#	info('main line')
#	p = Process(target=f, args=('bob',))
#	p.start()
#	p.join()

#import multiprocessing as mp
#import random
#import string

#output = mp.Queue()

#def rand_string(length, output):
#	rand_str = ''.join(random.choice( 
#			+ string.ascii_lowercase 
#			+ string.ascii_uppercase 
#			+ string.digits)
#			for i in range(length))
#	output.put(rand_str)

#processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]

#for p in processes:
#	p.start()
#

#for p in processes:
#	p.join()

#results = [output.get() for p in processes]

#print(results)

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
velvet_file_dict = test_annotate.annotate_velvet_contigs(args.RefSeq, args.output, args.directories, 20)












