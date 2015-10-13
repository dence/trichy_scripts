#Daniel Ence
#June 12, 2015

#open multiple tab-delimited taxonomer output files
#summarize read counts per TaxID
#assumes no ties (no special handling for ties)

from Bio import SeqIO
import re
import sys
#print 'Argument List:', str(sys.argv)

readsPerSeqID = []
taxIDs = {}
taxIDs.setdefault("unclassified",1)
#one file at a time
for file in sys.argv[2:len(sys.argv)]:	
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
for file in sys.argv[2:len(sys.argv)]:
	header = header + file + "\t"
header = header
#print str(taxIDs)

record_dict = SeqIO.to_dict(SeqIO.parse(open(sys.argv[1],'r'),"fasta"))


print header
for record in record_dict:
	taxID = record 
	line = taxID + "\t" + str(len(record_dict[taxID].seq)) + "\t"
	for currCounts in readsPerSeqID:
		if( currCounts.has_key(taxID)):
			#print currCounts
			line = line + str(currCounts[taxID]) + "\t"
		else:
			line = line + "0\t"
<<<<<<< HEAD
	line = line	
=======
>>>>>>> e695e7fedcdb804dd57736ca4506e5af4badc552
	print line
#print str(readsPerSeqID)

