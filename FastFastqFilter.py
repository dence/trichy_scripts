#Daniel Ence
#November 18, 2015

import abc
from Bio import SeqIO
from FastqFilter import FastqFilter
import re
import mmap


class FastFastqFilter(FastqFilter):

	def __init__(self):
		#self.oldFastq = fastqFile
		self.readsPerSeqID = {}
		self.readsFromOldFastq = {}
		self.newFastqList = []
		self.outputDirector = "."

	def main(fastqFileName,taxonomer,RefSeqs,outputDirectory):
		_init__()
		loadFastqFile(fastqFileName)
		loadTaxonomerOutput(taxonomer)
		filterFastq(RefSeqs,".")

	#fastfastq filter assumes that the fastq file has already been sorted. 
	def loadFastqFile(self, fastqFileName):
		fastq = open(fastqFilename,'r')
		
		fastqLengthMap = {}
		for readID,readBPs,plus,readQuals in fastq:
			print "readID is:\t" + readID + "\n"
				
		self.oldFastq = SeqIO.parse(fastqFileName,"fastq")

	def loadTaxonomerOutput(self, taxonomer):
		taxon_out = open(taxonomer,'r')
		for line in taxon_out:
			lineParts = re.split("\t",line,maxsplit=4)
			if(lineParts[0] == "C"):
				self.readsPerSeqID[lineParts[1]] = lineParts[3]
	
	def filterFastq(self, RefSeqs, outputDirectory):

		RefSeq_dict = SeqIO.to_dict(SeqIO.parse(open(RefSeqs,'r'),"fasta"))
		for record in self.oldFastq:
			if(record.id in self.readsPerSeqID and \
			self.readsPerSeqID[record.id] in RefSeq_dict):
				classifiedTaxa = self.readsPerSeqID[record.id]
				self.readsFromOldFastq.setdefault(classifiedTaxa,[]).append(record)

		for taxa in self.readsFromOldFastq:
			newFastq = outputDirectory + "/" + taxa + ".classified_reads.fastq"
			self.newFastqList.append(newFastq)
			newFastqHandle = open(newFastq,'w')
			currFastqList = self.readsFromOldFastq[taxa]
			SeqIO.write(currFastqList,newFastqHandle,"fastq")
			newFastqHandle.close()

	def getNewFastqList(self):
		return self.newFastqList	


if __name__ == "__main__":
	
