#Daniel Ence 
#08/10/2015


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import re
import subprocess


def gather_refseqs(refseq_fasta_file, sample_name):
	
	refseq_fasta_filehandle = SeqIO.parse(refseq_fasta_file,"fasta")
	sample_fasta_file = sample_name + "_velvet_contigs.fasta"
	sample_filehandle = open(sample_fasta_file,'a')

	for refseq in refseq_fasta_filehandle:
		contig_file = refseq.id + ".classified_reads.fastq_velvet/contigs.fa"	
		if(os.path.exists(contig_file)):
			contig_handle = SeqIO.parse(contig_file,"fasta")
			counter = 0
			for contig in contig_handle:
				curr_id = refseq.id + "_velvet_" + str(counter)
				tmp_new_record = SeqRecord(contig.seq,id=curr_id,name=curr_id,description="")
				SeqIO.write(tmp_new_record, sample_filehandle,"fasta")
				counter = counter + 1	
				

gather_refseqs(sys.argv[1], sys.argv[2])





