#Daniel Ence
#May 31, 2015
import sys
import re
import argparse


def main(RefSeqFile):
	dumb = FlatTaxonomyMaker()
	dumb.AddFlatTaxonomy(RefSeqFile)

class FlatTaxonomyMaker():
		
	def __init__(self):

	def AddFlatTaxonomy(self, FastaFileName):
		#writes a version of the fasta file with flat taxonomy
		#write a key file "key.txt" with the taxonomy too		

		fasta_file = open(FastaFileName,'r')
		
		key_file = open("key.txt",'w')
		added_fasta_file = "key_added_" + FastaFileName
		added_fasta = open(added_fasta_file,'w')

		fasta_pattern = re.compile("^\>")
		fasta_capture = re.compile("([\S]+)")
		taxID_count = 1
		root_string = "00__Root"
		root_line = str(taxID_count) + "\t0\t0\t" + root_string
		
		key_handle.write(root_line + "\n")
		for line in fasta_file:
			if(fasta_pattern.search(line)):
				match = fasta_capture.match(line,1)
				seq_ID = match.group(0)
				taxID_count = taxID_count + 1
				seq_ID_string = root_string + ";" + "{0:02d}".format(taxID_count) + "__" + seq_ID
				curr_line = str(taxID_count) + "\t1\t1\t" + seq_ID_string
				key_handle.write(curr_line + "\n")
				new_line = line.replace('\n',"\t" + seq_ID_string + "\n")
			else:
				fasta_hande.write(line)			
			
		key_handle.close()
		fasta_handle.close()
		fasta_file.close()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="adds a flat taxonomy to a fasta file")
	parser.add_argument("--RefSeq",type=str,help="Fasta file of sequences")
	args = parser.parse_args()
	main(args.RefSeq)


	
