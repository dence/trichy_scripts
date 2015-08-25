#Daniel Ence 
#08/10/2015


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import re
import subprocess
import argparse

import multiprocessing as mp

#take in file with the refseqs
#take in a list of the maker directories

#for each refseq
#gather the version of each refseq from each maker directory into a fasta file


def gather_est2genome_hits(refseq_fasta_file, gff_file):
	
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
				
def make_command_duples(RefSeq_file, maker_directories):

	print os.getcwd()	

	refseq_filehandle = SeqIO.parse(RefSeq_file, "fasta")
	
	duples = []
	for refseq in refseq_filehandle:
		curr_duple = (refseq, maker_directories)
		duples.append(curr_duple)
	return duples								

def find_longest_est2genome_hit(Seq_obj_name, curr_dir):

	datastore_file = re.sub("\.maker\.output","_datastore",curr_dir)
	curr_datastore = curr_dir + "/" + datastore_file 
	log_file = re.sub("\.maker\.output","_master_datastore_index.log",curr_dir)
	curr_datastore_log = curr_dir + "/" + log_file 
	
	new_record = SeqRecord("",id="",name="",description = "")
	longest_map = {}
	if(os.path.exists(curr_datastore_log)):
		log_handle = open(curr_datastore_log,'r')
		longest_hit = ""
		longest_length = 0
		for line in log_handle:
			seqname = line.split("\t")[0]
			if(re.search(Seq_obj_name,seqname) and re.search("FINISHED",line)):
				seq_dir = line.split("\t")[1]
				tmp_refseq = seq_dir.split("/")[len(seq_dir.split("/")) - 2].replace(".","%2E")#hardcoded this conversion
				gff_file = curr_dir + "/" + seq_dir + "/" +  tmp_refseq + ".gff"
				if(length_of_hit(gff_file) > longest_length):
					longest_map[gff_file] = length_of_hit(gff_file)
		log_handle.close()

		if(len(new_record.id) > 0 and has_est2genome_hit(longest_hit)):
			new_record = make_new_SeqRecord(longest_hit, curr_dir)

	return new_record

def has_est2genome_hit(gff_file):
	bool = 0
	gff_handle = open(gff_file, 'r')
	for gff_line in gff_handle:
		if(re.search("\texpressed_sequence_match\t",gff_line)):
			bool = 1
	return bool

def make_new_SeqRecord(gff_file, curr_sample_dir):
	hit_SeqID = ""
	hit_start = 0
	hit_stop = 0
	hit_strand = "+"
	gff_handle = open(gff_file, 'r')
	for gff_line in gff_handle:
		if(re.search("\texpressed_sequence_match\t",gff_line)):
			hit_SeqID = gff_line.split("\t")[0]
			hit_start = int(gff_line.split("\t")[3])
			hit_stop = int(gff_line.split("\t")[4])
			hit_strand = gff_line.split("\t")[6]
			
			if(hit_strand == "+"):
                                length = hit_stop - hit_start
                        else:
                                length = hit_start - hit_stop
                        break
        gff_handle.close()

	fasta_file = re.sub("\.maker\.output",".fasta",curr_sample_dir)
	fasta_handle = open(fasta_file,'r')
	record_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle,"fasta"))
	hit_record = record_dict[hit_SeqID]
	
	new_seq = hit_record.seq[hit_start - 1:hit_stop]
	if(hit_strand == "-"):
		new_seq = hit_record.seq[curr_start - 1:curr_stop].reverse_complement()
	
	seqname = curr_sample_dir.split("_")[0] + "_" + hit_SeqID
	return SeqRecord(new_seq,id=seqname,name=seqname,description="")

def length_of_hit(gff_file):
	#go through each line of the gff file
	#find the length of the est2genome hit
	
	gff_handle = open(gff_file, 'r')

	length = 0
	for gff_line in gff_handle:
		if(re.search("\texpressed_sequence_match\t",gff_line)):
			curr_start = int(gff_line.split("\t")[3])
			curr_stop = int(gff_line.split("\t")[4])
			curr_strand = gff_line.split("\t")[6]
		
			if(curr_strand == "+"):
				length = curr_stop - curr_start
			else:
				length = curr_start - curr_stop
			break
	gff_handle.close()

	return length 
	
def gather_and_align_est2genome_hits(command_duple):
	
	Seq_obj = command_duple[0]
	maker_directories = command_duple[1]		 				

	new_filehandle = open(Seq_obj.id + ".est2genome.fasta",'w')		
	
	SeqIO.write(Seq_obj,new_filehandle,"fasta")

	for curr_dir in maker_directories:
		#est2genome_Seq_obj = find_longest_est2genome_hit(Seq_obj.id, "../" + curr_dir)
		est2genome_Seq_obj = find_longest_est2genome_hit(Seq_obj.id, curr_dir)
		
		if(est2genome_Seq_obj.id != ""):
			SeqIO.write(est2genome_Seq_obj, new_filehandle,"fasta")

	new_filehandle.close()	

parser = argparse.ArgumentParser(description="gather est2genome hits from maker output directories")
parser.add_argument("RefSeq",type=str,help="fasta file of the reference sequences")
parser.add_argument("output_directory",help="directory to be made and filled with alignment files")
parser.add_argument("nprocs",type=int,help="number of processes to run concurrently")
parser.add_argument("directories",nargs='+',type=str,help="list of maker output directories")
args=parser.parse_args()

if not os.path.exists(args.output_directory):
	os.mkdir(args.output_directory)
#os.chdir(args.output_directory)

#command_duples = make_command_duples("../" + args.RefSeq, args.directories)
command_duples = make_command_duples(args.RefSeq, args.directories)

pool = mp.Pool(processes=args.nprocs)
map(gather_and_align_est2genome_hits, command_duples)
#pool.map(gather_and_align_est2genome_hits, command_duples)




