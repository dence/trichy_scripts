#Daniel Ence
#November 18, 2015

import abc
from Bio import SeqIO
from FastqFilter import FastqFilter
import re

class SlowFastqFilter(FastqFilter):

	def __init__(self):
		#self.oldFastq = fastqFile
		self.readsPerSeqID = {}
		self.readsFromOldFastq = {}
		self.newFastqList = []

	def loadFastqFile(self, fastqFileName):
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



