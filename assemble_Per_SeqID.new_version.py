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


def assemble_fastq(filteredFastq):

	velvetDir = filteredFastq + "_velvet"	
	velvetBin = "/home/dence/applications/velvet"
	velvetH_command = velvetBin + "/" + "velveth" + " " + velvetDir + " 29 " + " -fastq " + filteredFastq
	print velvetH_command
	subprocess.call(velvetH_command, shell=True)
	
	velvetG_command = velvetBin + "/" + "velvetg" + " " + velvetDir + " -cov_cutoff auto -exp_cov auto "
	print velvetG_command
	subprocess.call(velvetG_command,shell=True)

def prepare_RefSeq_dictionary(refseqFasta):
	RefSeqs_dictionary = {}
	RefSeqs = SeqIO.parse(refseqFasta, "fasta")
	for RefSeq in RefSeqs:
		RefSeqs_dictionary[RefSeq.id] = 1
	return RefSeqs_dictionary

def filter_fastqs(refseqs, taxonomerOut, fastq, prefix):
	print "refseqs is:\t" + refseqs
	print "taxonomerOut:\t" + taxonomerOut
	print "fastq is:\t" + fastq
	print "prefix is:\t" + prefix

	RefSeqs_dictionary = prepare_RefSeq_dictionary(refseqs)
	
	readsPerSeqID = {}
        taxon_out = open(taxonomerOut, 'r')
        for line in taxon_out:
		#print "line is:\t" + line
                lineParts = re.split("\t",line,maxsplit=4)
		#print line
		#print lineParts
                if(lineParts[0] == "C" and RefSeqs_dictionary.get(lineParts[3])):
                        re.compile("_")
                        IDParts = re.split("_", lineParts[1])
                        newID = str(" ").join(IDParts)
                        readsPerSeqID[IDParts[0]] = lineParts[3]
	
	print " reading from fastq: " + fastq + "\n"
	oldFastQ = SeqIO.parse(fastq, "fastq")
	readsFromOldFastq = {}
	recordCount = 0
	import time
	startFastqTime = time.clock()
	for record in oldFastQ:
		if(record.id in readsPerSeqID):
			currTaxa = readsPerSeqID[record.id]
			readsFromOldFastq.setdefault(currTaxa,[]).append(record)
	stopFastqTime = time.clock()
	elapse = stopFastqTime - startFastqTime
	print "Took this long to go through the fastq file:\t" + str(elapsed)
	

	newFastq_list = []
	for taxa in readsFromOldFastq:
		newFastq = args.prefix + "/" + taxa + ".classified_reads.fastq"
		newFastQ_list.append(newFastq)
		curr_handle = open(newFastq,'w')
		print "opened this fastq file:\t" + newFastq
		curr_fastq_list = readsFromOldFastq[taxa]
		print "Writing to this fastq file:\t" + newFastq
		SeqIO.write(curr_fastq_list, curr_handle, "fastq")
		curr_handle.close()

	return newFastq_list
	

def main():
	parser = argparse.ArgumentParser(description="Assemble reads classified to each SeqID")
	parser.add_argument("taxonomer_out",type=str, help="taxonomer output file")
	parser.add_argument("fastq", type=str, help="fastq file on which the taxonomer output is based")
	parser.add_argument("prefix", type=str, help="prefix for output files. A directory with this name will be created")
	parser.add_argument("RefSeqs", type=str, help="fasta file with reference sequences that we want to select reads and assemble on")
	parser.add_argument("steps",type=str,help="which steps to to run. Valid values are \"all\",\"fastq\",and\"assemble\"")
	parser.add_argument("nprocs", type=int, help="number of concurrent velvet jobs to run")
	args=parser.parse_args()

	print "taxonomer_out is:\t" + args.taxonomer_out
	print "fastq is:\t" + args.fastq
	print "prefix is:\t" + args.prefix
	print "RefSeqs is:\t" + args.RefSeqs
	print "steps is:\t" + args.steps
	print "nprocs is:\t" + str(args.nprocs)	
	
	fastq=0
	assemble=0
	if(args.steps == "all"):
		fastq=1
		assemble=1
	elif(args.steps == "fastq"):
	        fastq=1
	elif(args.steps == "assemble"):
		assemble=1
	else:
        	print "INVALID VALUE FOR STEPS!!!"
        	sys.exit()
	print "making directory for output " + args.prefix + "\n"
	#make a directory to hold the subsampled fastqs.
	if not os.path.exists(args.prefix):
	        os.mkdir(args.prefix)

	newFastQ_list = []
	if(fastq > 0):
		newFastQ_list = filter_fastqs(args.RefSeqs, args.taxonomer_out, args.fastq, args.prefix)
	else:
	#if we're not making fastqs, then get a list of what's already been made
		files = os.listdir(args.prefix)
		for file in files:
			if(file.endswith(".fastq")):
				newFastQ_list.append(args.prefix + "/" + file)

	if(assemble > 0):
		command_duples = []
		for fastQ in newFastQ_list:
			velvetDir = fastQ + "_velvet"
			mkdir_command = "mkdir " + velvetDir
			if not os.path.exists(velvetDir):
				print mkdir_command
				os.mkdir(velvetDir)
				velvet_to_run.append((fastQ))
			if not os.path.exists(velvetDir + "/" + "contigs.fa"):
				velvet_to_run.append((fastQ))

		pool = mp.Pool(processes=args.nprocs)
		#pool.map(assemble_fastq, command_duples) 
		map(assemble_fastq, command_duples) 

if __name__ == '__main__':
	main()

