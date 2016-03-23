#Daniel Ence 
#January 20, 2016

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
#import string
import os

parser = argparse.ArgumentParser(description="read in tabbed wublast output, output fasta sequences of the hit regions")
parser.add_argument("blastOut",type=str,help="tabbed wublast output")
parser.add_argument("fasta",type=str,help="fasta file used for blast query")
parser.add_argument("outFile",type=str,help="name of fasta file to create")
args=parser.parse_args()

#take in the blast output, map the seqID to the start and top of the blast hit 
#assumes only one blast hit per query/subject (B=1 V=1)

blastResults = open(args.blastOut,'r')

blastResultsMap = {}
for line in blastResults:
	lineParts = line.split("\t") 
	queryID = lineParts[0]
	queryStart = int(lineParts[17])
	queryEnd = int(lineParts[18])

	#print "query ID is:\t" + queryID
	#print "query start is:\t" + str(queryStart)
	#print "query end is:\t" + str(queryEnd)
	

	blastResultsMap[queryID] = (queryStart,queryEnd)

blastHits = []
fasta = SeqIO.parse(args.fasta,"fasta")
for query in fasta:

	hitSeq = ""
	hitName = ""

	#print "query is:\t" + query.id

	try:
		queryHitStart = blastResultsMap[query.id][0]
		queryHitEnd = blastResultsMap[query.id][1]
		hitSeq = query.seq[queryHitStart - 1:queryHitEnd]

		if(queryHitEnd < queryHitStart):
			#print "going here:"
			tmpSeq = query.seq[queryHitEnd:queryHitStart - 1]
			#print "tmpseq is:\t" + tmpSeq
			hitSeq = tmpSeq.reverse_complement()
			#print "hit seq is:\t" + hitSeq
		hitName = query.id
	except KeyError:
		hitSeq = query.seq
		hitName = query.id

	#print "query hit seq is:\t" + hitSeq			
		
	newRecord = SeqRecord(hitSeq, id=hitName,name=hitName,description="")
	blastHits.append(newRecord)
 	
outFastaFile = open(args.outFile,'w')
SeqIO.write(blastHits,outFastaFile,"fasta")
