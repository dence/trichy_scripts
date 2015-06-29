#Daniel Ence
#06/19/2015

import sys
import re


VCF = open(sys.argv[1], 'r')

FORMAT_found = 0
FORMAT_pos = 0
chrom_pattern = re.compile("CHROM")
hashtag_pattern = re.compile("^#")
FORMAT_pattern = re.compile("FORMAT")
header_line="overall_freq"
for line in VCF:
	line = line.rstrip('\n')
	if(FORMAT_pattern.search(line)):
		INFO_parts = re.split("\t", line)
		for part in INFO_parts:
			if(FORMAT_found == 1):
				header_line = header_line + "\t" + part	
			else:
				if FORMAT_pattern.match(part):
					FORMAT_found = 1
					FORMAT_pos = INFO_parts.index("FORMAT") + 1
	
	else:
		if(FORMAT_found == 1):
			varLineParts = re.split("\t", line)
			ref_allele = varLineParts[3]
			alt_allele = varLineParts[4]

			if(len(ref_allele) == 1 and len(alt_allele) == 1):
				curr_line = ""
				overall_INFO=varLineParts[FORMAT_pos - 2]
				overall_parts = re.split("\;", overall_INFO)
				overall_freq=overall_parts[1]
				overall_freq_parts=re.split("\=", overall_freq)
				overall_freq=overall_freq_parts[1]
				curr_line = curr_line + overall_freq + "\t"
				sample_INFOs = varLineParts[FORMAT_pos:len(varLineParts)]
				for sample_INFO in sample_INFOs:
					if(sample_INFO == "./."):
						curr_line = curr_line + "0\t"
					else:
						sample_parts = re.split("\:",sample_INFO)
						allele_depth = sample_parts[1]
						if(allele_depth != "."):
							[ref_depth,alt_depth] = re.split("\,",allele_depth)
							total=ref_depth + alt_depth
							alt_freq = float(alt_depth) / float(total)
							if(alt_freq > 0):
								curr_line = curr_line + "{0:.4f}".format(alt_freq) + "\t"
							else:
								curr_line = curr_line + "0\t"
						else:
							curr_line = curr_line + "0\t"
				print curr_line	
