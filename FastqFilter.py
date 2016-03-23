#Daniel Ence
#November 18, 2015

import abc

class FastqFilter(object):
	__metaclass__ = abc.ABCMeta

	#def __init__(self):
        #        oldFastq = ""
        #        readsPerSeqID = {}
        #        readsFromOldFastq = {}
        #        newFastqList = []

	@abc.abstractmethod
	def loadFastqFile(self, fastqFileName):
		"""Load a fastq file into memory"""
		return
	
	@abc.abstractmethod
	def loadTaxonomerOutput(self, taxonomer):
		"""Load a taxonomer file into memory"""
		return

	@abc.abstractmethod
	def filterFastq(self, RefSeqs, outputDirectory):
		"""Look for reads classified to RefSeqs and make a fastq for each RefSeq in /
		 the outputDirectory"""
		return

	def getNewFastqList(self):
		return self.newFastqList
