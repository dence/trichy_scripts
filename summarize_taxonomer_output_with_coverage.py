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

readsPerSeqID = []
taxIDs = {}
taxIDs.setdefault("unclassified",1)
#one file at a time
for file in args.taxonomer_results:	
	currFileCounts={}
	taxon_out = open(file, 'r')
	currTaxID = "unclassified"
	for line in taxon_out:
        	lineParts = re.split("\t",line,maxsplit=4)
        	if(lineParts[0] == "C"):
                	currTaxID=lineParts[3]
		currCount = currFileCounts.setdefault(currTaxID,0)
		currFileCounts[currTaxID] = currCount + 1
		taxIDs.setdefault(currTaxID,1)		
	readsPerSeqID.append(currFileCounts)
#make header line
header = "TaxID\tSeqLength\t"
for file in args.taxonomer_results: 
	header = header + file + "\t"
header = header
#print str(taxIDs)

record_dict = SeqIO.to_dict(SeqIO.parse(open(args.RefSeq,'r'),"fasta"))

counts_file = "taxonomer_counts.summary.txt"
coverage_file = "taxonomer_coverage.summary.txt"
summary_counts = open(counts_file,'w')
summary_coverage = open(coverage_file,'w')
summary_counts.write(header + "\n")
summary_coverage.write(header + "\n")
#print header
for record in record_dict:
	taxID = record 
	line = taxID + "\t" + str(len(record_dict[taxID].seq)) + "\t"
	cov_line = taxID + "\t" + str(len(record_dict[taxID].seq)) + "\t"
	
	for currCounts in readsPerSeqID:
		if( currCounts.has_key(taxID)):
			#print currCounts
			line = line + str(currCounts[taxID]) + "\t"
			#assuming 100 bp reads here.
			cov = float( float(currCounts[taxID]) * float(100)) / float(len(record_dict[taxID].seq))
			cov_string = '{:.6}'.format(cov)
			cov_line = cov_line + cov_string + "\t"
		else:
			line = line + "0\t"
			cov_line = cov_line + "0\t"
	#print line
	summary_counts.write(line + "\n")
	summary_coverage.write(cov_line + "\n")
#print str(readsPerSeqID)

