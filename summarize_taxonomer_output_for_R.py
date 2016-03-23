#Daniel Ence
#June 12, 2015

#open multiple tab-delimited taxonomer output files
#summarize read counts per TaxID
#assumes no ties (no special handling for ties)

from Bio import SeqIO
import argparse
import re
import sys
#print 'Argument List:', str(sys.argv)

parser = argparse.ArgumentParser(description="Summary read count and effective coverage of a set of taxonomer results")
parser.add_argument("RefSeq",type=str,help="File of Sequences to get results for")
parser.add_argument("taxonomer_results",nargs='+',type=str,help="list of taxonomer outputs (assumes no ties)")
args = parser.parse_args()

RefSeqDict = SeqIO.to_dict(SeqIO.parse(open(args.RefSeq,'r'),"fasta"))

#print out a report where each each row is an effective coverage measurment
#with the columns: sample name, gene name, effective coverage for that gene in that sample

header = "Gene Name\tSample Name\tEffective Coverage\n"
coverage_file = "taxonomer_coverage_for_R.summary.txt"
summary_coverage = open(coverage_file,'w')
summary_coverage.write(header)

for file in args.taxonomer_results:
	currFileCounts = {}
	currTaxonomerOutput = open(file,'r')

	taxIDs = {}
		
	currentGeneName = "unclassified"
	for line in currTaxonomerOutput:
		lineParts = re.split("\t",line,maxsplit=4)
		if(lineParts[0] == "C"):
			currentGeneName = lineParts[3]
		currCount = currFileCounts.setdefault(currentGeneName,0)
		currFileCounts[currentGeneName] = currCount + 1
		taxIDs.setdefault(currentGeneName,1)

	for record in RefSeqDict:
		
		cov_line = record + "\t" + file + "\t"
		cov_string = "0"
		if(currFileCounts.has_key(record)):
			cov = float( float(currFileCounts[record]) * float(100)) / float(len(RefSeqDict[record].seq))
			
			#if(cov < 1):
			#	cov_string = '{:.6}'.format(cov)
			#else:
			#	cov_string = str(cov)
			
			cov_string = str(cov)
			
		cov_line = cov_line + cov_string + "\n"	
		summary_coverage.write(cov_line)
			

