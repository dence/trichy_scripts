#Daniel Ence
#June 12, 2015

#open multiple tab-delimited taxonomer output files
#summarize read counts per TaxID
#assumes no ties (no special handling for ties)

import re
import sys
#print 'Argument List:', str(sys.argv)

readsPerSeqID = []
taxIDs = {}
taxIDs.setdefault("unclassified",1)
#one file at a time
for file in sys.argv[1:len(sys.argv)]:	
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
header = "TaxID\t"
for file in sys.argv[1:len(sys.argv)]:
	header = header + file + "\t"
header = header + "\n"
#print str(taxIDs)
print header
for taxID in taxIDs:
	line = taxID + "\t"
	for currCounts in readsPerSeqID:
		if( currCounts.has_key(taxID)):
			#print currCounts
			line = line + str(currCounts[taxID]) + "\t"
		else:
			line = line + "0\t"
	line = line + "\n"	
	print line
#print str(readsPerSeqID)

