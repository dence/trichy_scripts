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


def gather_cdna2genome_hits(refseq_fasta_file, gff_file):
	
	refseq_fasta_filehandle = SeqIO.parse(refseq_fasta_file,"fasta")
	sample_fasta_file = sample_name + "_velvet_contigs.fa"
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

def find_longest_cdna2genome_hit(Seq_obj_name,curr_sample_gff):
	gff_handle = open(curr_sample_gff,'r')

	new_record = SeqRecord("",id="",name="",description="")
	longest_hit = ""
	longest_length = 0
	for line in gff_handle:
		seqname = line.split("\t")[0]
		if(re.search(Seq_obj_name,seqname) and re.search("\texpressed_sequence_match\t",line)):
			print "matched!!"
			print line
			curr_start = int(line.split("\t")[3])
                        curr_stop = int(line.split("\t")[4])
			length_of_hit = curr_start - curr_stop
			if(length_of_hit < 0):
				length_of_hit = length_of_hit * -1
			if(length_of_hit > longest_length):
				longest_length = length_of_hit
				longest_hit = line
	if(longest_length > 0):
			print "longest_cdna2genome hit for " + Seq_obj_name + " in this file:\t" + curr_sample_gff
			print longest_hit
			curr_sample_fasta = curr_sample_gff.replace("all.gff","fasta")
			new_record = make_new_SeqRecord(longest_hit,curr_sample_fasta)
	else:
		print "no cdna2genome hit for " + Seq_obj_name + " in this file:\t" + curr_sample_gff 
	return new_record

def make_new_SeqRecord(gff_line, curr_sample_fasta):
	print "making a SeqRecord from this line:\t"
	print gff_line
	hit_SeqID = gff_line.split("\t")[0]
	hit_start = int(gff_line.split("\t")[3])
	hit_stop = int(gff_line.split("\t")[4])
	hit_strand = gff_line.split("\t")[6]
			
	if(hit_strand == "+"):
		length = hit_stop - hit_start
	else:
		length = hit_start - hit_stop

	fasta_handle = open(curr_sample_fasta,'r')
	hit_record = ""
	print "opening this file:\t" + curr_sample_fasta

	#fasta_dict = SeqIO.to_dict(SeqIO.parse(open(curr_sample_fasta,'r'),"fasta"))
	#hit_record = fasta_dict[hit_SeqID]
	for record in SeqIO.parse(fasta_handle,"fasta"):
		if(record.id == hit_SeqID):
			hit_record = record
	
	new_seq = hit_record.seq[hit_start - 1:hit_stop]
	if(hit_strand == "-"):
		new_seq = hit_record.seq[hit_start - 1:hit_stop].reverse_complement()
	
	seqname = curr_sample_fasta.split("_")[0] + "_" + hit_SeqID
	return SeqRecord(new_seq,id=seqname,name=seqname,description="")

def length_of_hit(gff_file):
	#go through each line of the gff file
	#find the length of the cdna2genome hit
	
	gff_handle = open(gff_file, 'r')

	length = 0
	for gff_line in gff_handle:
		if(re.search("\texpressed_sequence_match\t",gff_line)):
			curr_start = int(gff_line.split("\t")[3])
			curr_stop = int(gff_line.split("\t")[4])
			curr_strand = gff_line.split("\t")[6]
			length = curr_stop - curr_start	
			break
	gff_handle.close()

	if(length < 0):
		length = length * -1
	return length 
	
def gather_and_align_cdna2genome_hits(command_duple):
	
	Seq_obj = command_duple[0]
	maker_directories = command_duple[1]		 				

	cdna2genome_fasta = Seq_obj.id + ".cdna2genome.fasta"
	new_filehandle = open(cdna2genome_fasta,'w')		
	print "writing to this file:\t" + cdna2genome_fasta
	
	SeqIO.write(Seq_obj,new_filehandle,"fasta")

	for curr_gff in maker_directories:
		cdna2genome_Seq_obj = find_longest_cdna2genome_hit(Seq_obj.id,curr_gff)

		if(cdna2genome_Seq_obj.id != ""):
			SeqIO.write(cdna2genome_Seq_obj,new_filehandle,"fasta")

	new_filehandle.close()
	print "finished with this file:\t" + cdna2genome_fasta

	#now align the new cdna2genome hits file
	curr_align_out = Seq_obj.id + ".cdna2genome.aln"
	clustalw2_string = "clustalw2 -INFILE=" + cdna2genome_fasta + " -align -type=DNA -outfile=" + curr_align_out
	print "aligning with this command:\t" + clustalw2_string	
	subprocess.call(clustalw2_string,shell=True)
	print "finished_aligning this:\t" + curr_align_out

parser = argparse.ArgumentParser(description="gather cdna2genome hits from maker output directories")
parser.add_argument("RefSeq",type=str,help="fasta file of the reference sequences")
#parser.add_argument("output_directory",help="directory to be made and filled with alignment files")
parser.add_argument("nprocs",type=int,help="number of processes to run concurrently")
parser.add_argument("directories",nargs='+',type=str,help="list of maker output directories")
args=parser.parse_args()

#if not os.path.exists(args.output_directory):
#	os.mkdir(args.output_directory)
#os.chdir(args.output_directory)

#command_duples = make_command_duples("../" + args.RefSeq, args.directories)
command_duples = make_command_duples(args.RefSeq, args.directories)

#pool = mp.Pool(processes=args.nprocs)
map(gather_and_align_cdna2genome_hits, command_duples)
#pool.map(gather_and_align_cdna2genome_hits, command_duples)




