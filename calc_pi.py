#Daniel Ence
#06/19/2015

import sys
import re


VCF = open(sys.argv[1], 'r')

pi_counts = {}
sample_positions = []
FORMAT_found = 0
FORMAT_pos = 0
chrom_pattern = re.compile("CHROM")
hashtag_pattern = re.compile("^#")
FORMAT_pattern = re.compile("FORMAT")
for line in VCF:
	line = line.rstrip('\n')
	if(FORMAT_pattern.search(line)):
		INFO_parts = re.split("\t", line)
		for part in INFO_parts:
			if(FORMAT_found == 1):
				sample_positions.append(part)
				pi_counts.setdefault(part, [0,0])
			else:
				if FORMAT_pattern.match(part):
					FORMAT_found = 1
					FORMAT_pos = INFO_parts.index("FORMAT") + 1
	else:
		if(FORMAT_found == 1):
			varLineParts = re.split("\t", line)
			for sample in pi_counts:
				sample_pos = FORMAT_pos + sample_positions.index(sample)
				sample_INFO = varLineParts[sample_pos]
				ref_allele = varLineParts[3]
				alt_allele = varLineParts[4]
				if(sample_INFO != "./." and
				(len(ref_allele) == 1 and len(alt_allele) == 1)): 
					INFO_parts = re.split("\:", sample_INFO)
					curr_GT = INFO_parts[0]
					#print curr_GT
					if(curr_GT == "0/0"):
						GT_list = pi_counts[sample]
						GT_list[1] = GT_list[1] + 1
					else:
						if(curr_GT == "0/1"):
							GT_list = pi_counts[sample]
							GT_list[1] = GT_list[1] + 1
							GT_list[0] = GT_list[0] + 1
						else:
							GT_list = pi_counts[sample]
							GT_list[1] = GT_list[1] + 1
							GT_list[0] = GT_list[0] + 1
					#print pi_counts
for sample in pi_counts:
	curr_pi_list = pi_counts[sample]
	diffs = curr_pi_list[0]
	positions = curr_pi_list[1]	
	pi = float(curr_pi_list[0]) / float(curr_pi_list[1])
	print sample + "\t" + str(diffs) + "\t" + str(positions) + "\t" + "{0:.4f}".format(pi)
