

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

args=parser.parse_args()

print "making directory for output " + args.prefix + "\n"
#make a directory to hold the subsampled fastqs.
if not os.path.exists(args.prefix):
	os.mkdir(args.prefix)


print "reading from " + args.taxonomer_out + "\n"
readsPerSeqID = {}
SeqIDsFileHandles = {}
newFastQ_list = []
taxon_out = open(args.taxonomer_out, 'r')
for line in taxon_out:
	lineParts = re.split("\t",line,maxsplit=4)
	if(lineParts[0] == "C"):
		re.compile("_")
		IDParts = re.split("_",lineParts[1])
		newID = str(" ").join(IDParts)
		print "newID is:\t" + newID
		readsPerSeqID[IDParts[0]] = lineParts[3]
		if(lineParts[3] in SeqIDsFileHandles):
			#ddooo
			tmp = 1	
		else:
			handleName = args.prefix + "/" + lineParts[3] + ".classified_reads.fastq"
			newFastQ_list.append(handleName)
			curr_handle = open(handleName,"w")
			print "opened this file handle:\t" + handleName + "\n"	
			SeqIDsFileHandles[lineParts[3]] = curr_handle

#foreach key in reads_per_seqID
#open the old fastqs, subsample based on the taxonomer output
#
print "	reading from fastq: " + args.fastq + "\n"
oldFastQ = SeqIO.parse(args.fastq, "fastq")
recordCount = 0
for record in oldFastQ:
	recordCount = recordCount + 1
	print "read this many records:\t" + str(recordCount)
	print "recor.id is:\t" + record.id + "\n"
	if(record.id in readsPerSeqID):
		print "matched to " + record.id
		curr_handle = SeqIDsFileHandles[readsPerSeqID[record.id]]
		SeqIO.write([record], curr_handle, "fastq")	
			
for seqID in SeqIDsFileHandles:
	curr_handle = SeqIDsFileHandles[seqID]
	curr_handle.close()
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
	velvetH_command = velvetBin + "/" + "velveth" + " " + velvetDir + " 29 " + " -fastq " + fastQ 
	print velvetH_command
	subprocess.call(velvetH_command,shell=True)
	velvetG_command = velvetBin + "/" + "velvetg" + " " + velvetDir + " -cov_cutoff auto -exp_cov auto "
	print velvetG_command
	subprocess.call(velvetG_command,shell=True)
print "finished\n"	
