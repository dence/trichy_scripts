#Daniel Ence 
#06/11/2015

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

#read in the taxonomer output
#put read IDs in a hash, mapped to the seqID to which they were classified

parser = argparse.ArgumentParser(description="Assemble reads classified to each SeqID")
parser.add_argument("taxonomer_out",type=str, help="taxonomer output file")
parser.add_argument("fastq", type=str, help="fastq file on which the taxonomer output is based")
parser.add_argument("prefix", type=str, help="prefix for output files. A directory with this name will be created")
parser.add_argument("RefSeqs", type=str, help="fasta file with reference sequences that we want to select reads and assemble on")
parser.add_argument("steps",type=str,help="which steps to to run. Valid values are \"all\",\"fastq\",and\"assemble\"")
args=parser.parse_args()

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
	RefSeqs_dictionary = {}
	RefSeqs = SeqIO.parse(args.RefSeqs,"fasta")
	for RefSeq in RefSeqs:
		RefSeqs_dictionary[RefSeq.id] = 1

	print "reading from " + args.taxonomer_out + "\n"
	readsPerSeqID = {}
	#SeqIDsFileHandles = {}
	#newFastQ_list = []
	taxon_out = open(args.taxonomer_out, 'r')
	for line in taxon_out:
		lineParts = re.split("\t",line,maxsplit=4)
		if(lineParts[0] == "C" and RefSeqs_dictionary.get(lineParts[3])):
			re.compile("_")
			IDParts = re.split("_",lineParts[1])
			newID = str(" " ).join(IDParts)
			#print "this read:\t" + IDParts[0] + "\t to this taxID\t" + lineParts[3]
			readsPerSeqID[IDParts[0]] = lineParts[3]

	#foreach key in reads_per_seqID
	#open the old fastqs, subsample based on the taxonomer output
	#
	print "	reading from fastq: " + args.fastq + "\n"
	oldFastQ = SeqIO.parse(args.fastq, "fastq")
	readsFromOldFastq = {}
	recordCount = 0
	import time
	startFastqTime = time.clock()
	for record in oldFastQ:
	#	print "read this many records:\t" + str(recordCount)
	#	print "record.id is:\t" + record.id + "\n"
		if(record.id in readsPerSeqID):
	#		print "matched to this ID:\t" + record.id
			currTaxa = readsPerSeqID[record.id]
			readsFromOldFastq.setdefault(currTaxa,[]).append(record)
		#recordCount = recordCount + 1
	stopFastqTime = time.clock()
	elapsed = stopFastqTime - startFastqTime
	print "Took this long to go through the fastq file:\t" + str(elapsed)
	newFastQ_list = []
	for taxa in readsFromOldFastq:
		new_fastq = args.prefix + "/" + taxa + ".classified_reads.fastq"
		newFastQ_list.append(new_fastq)
		curr_handle = open(new_fastq,'w')
		print "opened this fastq file:\t" + new_fastq
		curr_fastq_list = readsFromOldFastq[taxa]
		print "Writing to this fastq file:\t" + new_fastq
		SeqIO.write(curr_fastq_list,curr_handle,"fastq")
		curr_handle.close()
		print "Finished writing to this fastq file:\t" + new_fastq
else:
	#if we're not making new fastqs, get the list of what's already been  made
	files = os.listdir(args.prefix)
	for file in files:
		if(file.endswith(".fastq")):
			newFastQ_list.append(args.prefix + "/" + file)
#assemble the new fastq with velvet
#hardcode a path to velvet here
for fastQ in newFastQ_list:
	print "running velvet for " + fastQ + "\n"
	velvetDir=fastQ + "_velvet"
	velvetBin = "/home/dence/applications/velvet/"
	mkdir_command = "mkdir " + velvetDir
	if not os.path.exists(velvetDir):
		print mkdir_command
	        os.mkdir(velvetDir)
	#check whether velvet finished running
	if not os.path.exists(velvetDir + "/" + "contigs.fa"):
		velvetH_command = velvetBin + "/" + "velveth" + " " + velvetDir + " 29 " + " -fastq " + fastQ 
		print velvetH_command
		subprocess.call(velvetH_command,shell=True)
		velvetG_command = velvetBin + "/" + "velvetg" + " " + velvetDir + " -cov_cutoff auto -exp_cov auto "
		print velvetG_command
		subprocess.call(velvetG_command,shell=True)
	else:
		print "Finished running velvet for:\t" + velvetDir + "/" + "contigs.fa"

			

print "finished\n"
	
