#Daniel Ence
#06/20/2015

import sys
import re
import os
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def     get_longest_contig(fasta_filename):
        curr_longest_contig = None
        len_curr_longest_contig = -1
        for contig in SeqIO.parse(fasta_filename, "fasta"):     
                if(len(contig.seq) > len_curr_longest_contig):
                        curr_longest_contig = contig
                        len_curr_longest_contig = len(contig)
        return curr_longest_contig

#take a fasta file of reference sequences
#and a list from the commandline of directories
#that contain velvet assemblies of reads that classified
#to the reference sequences
refSeq_file = sys.argv[1]
print "File with Referenc sequences is:\t" + refSeq_file
assembly_directories = sys.argv[2:len(sys.argv)]
print "assembly directories is:\t"
print  assembly_directories

#make a new directory for the alignments
#for now call directory "test_alignment"
if not os.path.exists("test_alignment"):
	print "making this directory:\ttest_alignment\n"
	os.mkdir("test_alignment")

#make a fasta file with the refseq and the 
#velvet results for the refSeq
refSeqs = SeqIO.parse(refSeq_file, "fasta")
for refSeq in refSeqs:
	curr_SeqID = refSeq.id
	print "curr_SeqID\t" + curr_SeqID						
	
	curr_SeqID_fasta = "test_alignment" + "/" + curr_SeqID + "_alignment.fasta" 
	curr_fasta_handle = open(curr_SeqID_fasta, 'w')
	
	SeqIO.write([refSeq], curr_fasta_handle, "fasta")	
	
	for dir in assembly_directories:	
		curr_velvet_fasta = dir + "/" + refSeq.id + ".classified_reads.fastq_velvet" + "/" + "contigs.fa"
		print "curr velvet contig is:\t" + curr_velvet_fasta

		if(os.path.isfile(curr_velvet_fasta) and \
		os.stat(curr_velvet_fasta).st_size > 0):	
			#read first entry (should only be one)
			print "reading from thsi file:\t" + curr_velvet_fasta

			curr_velvet_contig = get_longest_contig(curr_velvet_fasta)
			#curr_velvet_contig = SeqIO.read(open(curr_velvet_fasta,'r'), "fasta")
			new_ID = dir[2:len(dir)] + "_velvet"
			new_record = SeqRecord(Seq(str(curr_velvet_contig.seq)),id=new_ID,name=curr_velvet_contig.name, \
			description=curr_velvet_contig.description)
			SeqIO.write(new_record, curr_fasta_handle, "fasta")
	
	curr_fasta_handle.close()
	#align with clustalw2
	curr_seqID_align_out = "test_alignment" + "/" + curr_SeqID + ".pp.txt"
	clustalw2_command = "clustalw2 -INFILE=" + curr_SeqID_fasta + " -align -type=DNA -outfile=" + curr_seqID_align_out \
	+ " -output=phylip-relaxed"
	print "aliging this file:\t" + curr_SeqID_fasta
	print clustalw2_command
	#subprocess.call(clustalw2_command,shell=True) 


					
