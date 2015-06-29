#Daniel Ence

import sys
import re
fasta_file = open(sys.argv[1], 'r')

key_handle = open("key.txt",'w')
fasta_handle = open("key_added_" + sys.argv[1], 'w')

fasta_pattern = re.compile("^\>")
fasta_capture = re.compile("([\S]+)")
taxID_count = 1
root_string = "00__Root"
root_line = str(taxID_count) + "\t0\t0\t" + root_string
key_handle.write(root_line)
for line in fasta_file:
	if(fasta_pattern.search(line)):
		match = fasta_capture.match(line,1)
		seq_ID = match.group(0)
		taxID_count = taxID_count + 1
		seq_ID_string = root_string + ";" + "{0:02d}".format(taxID_count) + "__" + seq_ID
		curr_line = str(taxID_count) + "\t1\t1\t" + seq_ID_string
		key_handle.write(curr_line)
		new_line = line.replace('\n',"\t" + seq_ID_string + "\n")
		fasta_handle.write(new_line)
	else:
		fasta_handle.write(line)
			
key_handle.close()
fasta_handle.close()
fasta_file.close()

		
