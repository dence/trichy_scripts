#Daniel Ence 
#01/25/2016

#read in the taxonomer output
#put read IDs in a hash, mapped to the seqID to which they were classified

#select those reads from a fastq (use external script??)
#output to new fastq

#Assemble the new fastq with velvet (velveth, then velvetg)

#arguments to take in:
#taxonomer output
#fastqfile
#prefix for outputs/intermediate files
#

import argparse
import sys
import os
import subprocess
from Bio import SeqIO
import re

import multiprocessing as mp

import glob

import pdb

def main(outputDirectory,cpus):
	tmp = assemble_fastqs(outputDirectory)
	tmp.assemble_batch_fastqs(cpus)

class assemble_fastqs():

	def __init__(self, outputDirectory):
		self.outputDirectory = outputDirectory

	def assemble_batch_fastqs(self,cpus):
		print "here outputDirectory is:\t" + self.outputDirectory
		print "here cpus is:\t" + cpus
		#get a list of files in the output directory
		fastqs = glob.glob( self.outputDirectory + "/*classified_reads.fastq")
		#fastqs = glob.glob("outputDirectory" + "/")
		#print "fastqs is:\t"
		#print fastqs

		#make a list of command duples
		#run N(=cpus) number of assembly procs
                pool = mp.Pool(processes=int(cpus))
                pool.map(self.assemble_fastq, fastqs)
                #map(self.assemble_fastq, fastqs)
	
	def assemble_fastq(self, fastq):
		print "running velvet for " + fastq + "\n"
		velvetDir = fastq + "_velvet"
		velvetBin = "/home/dence/applications/velvet/"
		mkdir_command = "mkdir " + velvetDir

		if not os.path.exists(velvetDir):
			print mkdir_command
			os.mkdir(velvetDir)
		
		if not os.path.exists(velvetDir + "/" + "contigs.fa"):
			velvetH_command = velvetBin + "/" + "velveth" + " " + velvetDir + " 29 " + " -fastq " + fastq
			print velvetH_command
			subprocess.call(velvetH_command, shell=True)
			velvetG_command = velvetBin + "/" + "velvetg" + " " + velvetDir + " -cov_cutoff auto -exp_cov auto "
			print velvetG_command
			subprocess.call(velvetG_command,shell=True)
		else:
			"Finished running velvet for:\t" + velvetDir + "/" + "contigs.fa"		

if __name__ == "__main__":
	main(sys.argv[1],sys.argv[2])
