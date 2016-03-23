#Daniel Ence
#November 18, 2015

import abc
from Bio import SeqIO
from FastqFilter import FastqFilter
import re
import mmap
import sys
import os

def main(fastqFileName,taxonomer,RefSeqs,outputDirectory):
	tmp = FastFastqFilter()
	tmp.loadFastqFile(fastqFileName)
	tmp.loadTaxonomerOutput(taxonomer,RefSeqs)
	tmp.filterFastq(RefSeqs,outputDirectory)


class FastFastqFilter(FastqFilter):

	def __init__(self):
		self.oldFastq = ""
		self.readsPerSeqID = {}
		self.readsFromOldFastq = {}
		self.newFastqList = []
		self.ReadMap = {}

	#fastfastq filter assumes that the fastq file has already been sorted. 
	def loadFastqFile(self, fastqFileName):

		self.oldFastq = fastqFileName
		byte_start = 0
		N=4
		fastqFile = open(fastqFileName,'r')
		lines = []
		for line in fastqFile:
			lines.append(line)
			if len(lines) == N:
				byte_length = 0
				byte_length = len(lines[0])
				#print "in loadFastqFile, readID is:\t" + lines[0]
				byte_length = byte_length + len(lines[1])
				byte_length = byte_length + len(lines[2])
				byte_length = byte_length + len(lines[3])
				match = re.search('\@(\S+)',lines[0])
				readID = match.group()[1:]
				#print "in loadFastqFile, readID is:\t" + readID
				self.ReadMap[readID] = (byte_start,byte_length)
				byte_start = byte_start + byte_length
				lines = []
					
		fastqFile.close()
		#print self.byteMap 	
		#store the read ID 
		#map it to a duple with 
		#the byte start and byte length for this entry

	def loadTaxonomerOutput(self, taxonomer,RefSeqs):
		
		RefSeqs_dictionary = {}
		RefSeqs = SeqIO.parse(RefSeqs,"fasta")
		for RefSeq in RefSeqs:
			RefSeqs_dictionary[RefSeq.id] = 1

		taxon_out = open(taxonomer,'r')
		for line in taxon_out:
			lineParts = re.split("\t",line,maxsplit=4)
			if(lineParts[0] == "C" and RefSeqs_dictionary.get(lineParts[3])):
				IDParts = re.split("_",lineParts[1])
				newID = str(" ").join(IDParts)
				self.readsPerSeqID[lineParts[1]] = lineParts[3]
		sortedReadsPerSeqID = [ (k,v) for v,k in sorted([(v,k) for k,v in self.readsPerSeqID.items()])]
		self.readsPerSeqID = sortedReadsPerSeqID	
	
	def filterFastq(self, RefSeqs, outputDirectory):
	
		os.mkdir(outputDirectory)	
		for readID,RefSeq in self.readsPerSeqID:
			#print "in filter fastq readID is:\t" + readID
			start,length = self.ReadMap[readID]
			
			fastqFile = open(self.oldFastq,'r')
			tmpMap = mmap.mmap(fastqFile.fileno(),0,access=mmap.ACCESS_READ)
			tmpMap.seek(start)
			fastqEntry = tmpMap.read(length)
			fastqFile.close()

			newFastq = outputDirectory + "/" + RefSeq + ".classified_reads.fastq"
			self.newFastqList.append(newFastq)# need to make this a map and prevent dups
			newFastqHandle = open(newFastq,'a')
			newFastqHandle.write(fastqEntry)
			newFastqHandle.close()	

	def getNewFastqList(self):
		return self.newFastqList	

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])	
