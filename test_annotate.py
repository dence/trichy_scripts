from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import re
import subprocess

import multiprocessing as mp


def write_contigs(to_filehandle,from_filehandle,curr_sample,curr_refseq_id):
	
	counter = 0
	for contig in from_filehandle:
		curr_id = curr_sample + "_" + str(counter) + "." + curr_refseq_id
		counter = counter + 1
		tmp_record = SeqRecord(contig.seq,id=curr_id,name=curr_id,description="")
		SeqIO.write(tmp_record,to_filehandle,"fasta")

def gather_refseq_contigs(refseq_obj,output_dir,directories):
	
	new_fasta_file = output_dir + "/" + refseq_obj.id + ".velvet_contigs.fa"
	new_filehandle = open(new_fasta_file, 'w')

	for curr_dir in directories:
		curr_velvet_file = curr_dir + "/" + refseq_obj.id + ".classified_reads.fastq_velvet/contigs.fa"
		curr_sample = re.split("/",curr_dir)[1]	
		if(os.path.exists(curr_velvet_file)):
			
			#can't use the dictionary method because 
			#we sometimes get multiple contigs
			#curr_fasta_dict = SeqIO.to_dict(SeqIO.parse(curr_velvet_file,"fasta"))
			curr_fasta_file = SeqIO.parse(curr_velvet_file,"fasta")
			write_contigs(new_filehandle,curr_fasta_file,curr_sample,refseq_obj.id)
		else:
			print "This file was not found:\t" + curr_velvet_file + " in this directory:\t" + curr_dir
	new_filehandle.close()

	return new_fasta_file

def run_maker(filename):
	
	maker_command_string = "maker -g " + filename + " " + tmp_opts + " maker_bopots.ctl maker_exe.ctl"
	subprocess(maker_string,shell=True)

def annotate_velvet_contigs(RefSeq_file,output_dir,directories,number_procs):
	
	#figuring out which files go with which refseqs	
	refseq_fasta = SeqIO.parse(RefSeq_file,"fasta")
	
	pool = mp.Pool(processes=number_procs)

	results = [pool.apply( gather_refseq_contigs, args=(refseq,output_dir,directories)) for refseq in refseq_fasta]
	refseq_fasta.close()

	subprocess.call("cd " + output_dir,shell=True)
	subprocess.call("maker -CTL",shell=True)
	
	opts_handle = open("maker_opts.ctl",'r')
	tmp_opts = "tmp_opts"
	tmp_opts_handle = open("tmp_opts",'w')

	for line in opts_handle:
		if(re.match("est=",line):
			tmp_opts_handle.write(line.replace("est= ","est=" + RefSeq_file + " "))
		elif re.match("single_exon=",line):
			tmp_opts_handle.write(line.replace("single_exon=0","single_exon=1"))
		elif re.match("est2genome=",line):
			tmp_opts_handle.write(line.replace("est2genome=0","est2genome=1"))
		elif re.match("model_org=",line):
			tmp_opts_handle.write(line.replace("model_org=all","model_org="))
		elif re.match("always_complete=",line):
			tmp_opts_handle.write(line.replace("always_complete=0","always_complete=1"))
		else:
			tmp_opts_handle.write(line)
	opts_handle.close()
	tmp_opts_handle.close()

	results = [pool.appy( run_maker, args=(new_file)) for new_file in results]
	
def gather_est2genome_seqs(refseq_obj, est2genome_handle, log_line, velvet_file):
	seq_dir = log_line.split("\t")[1]
	tmp_refseq = seq_dir.split("/")[3].replace(".","%2E")#hardcoded in this position
	gff_file = refseq_obj.id + ".velvet_contigs.maker.output/" + seq_dir + "/" + tmp_refseq + ".gff"
	gff_handle = open(gff_file,'r')
	for gff_line in gff_handle:
		if(re.search("est2gneome",gff_line) and \
		re.search("\texpressed_sequence_match\t",gff_line)):
			curr_start = int(gff_line.split("\t")[3])
			curr_stop = int(gff_line.split("\t")[4])
			curr_strand = gff_line.split("\t")[6]
			
			tmp_handle = open(velvet_file,'r')
			tmp_fasta = SeqIO.to_dict(SeqIO.parse(tmp_handle,"fasta"))
			tmp_handle.close()
		
			if seq_dir.split("/")[3] in tmp_fasta:
				curr_record = tmp_fasta[seq_dir.split("/")[3]]
			else:
				continue
			new_seq = curr_record.seq[curr_start - 1:curr_stop]
			if(curr_strand == "-"):
				new_seq = curr_record.seq[curr_start - 1:curr_stop].reverse_complement()
			new_record = SeqRecord(new_seq,id=seqname,name=seqname,description="")
				
			SeqIO.write(est2genome_handle,"fasta")
	
def gather_est2genome_hits(refseq_obj, output_dir, velvet_file):
	
	est2genome_fasta = output_dir + "/" + refseq.id + ".est2genome.fasta"
	est2genome_handle = open(est2genome_fasta,'w')

	curr_datastore = refseq_obj.id + ".velvet_contigs.maker.output/" + refseq_obj.id \
	+ ".velvet_contigs_datastore"
	curr_datastore_log = refseq.id + ".velvet_contigs.maker.output/" + refseq_obj.id \
	+ ".velet_contigs_master_datastore_index.log"
	
	log_handle = open(curr_datastore_log,'r')
	for line in log_handle:
		seqname = line.split("\t")[0]
		if(re.search(refseq.id,seqname) and re.search("FINISHED",line)):
			gather_est2genome_seqs(refseq_obj, est2genome_handle,line, velvet_file)
	log_handle.close()	

def run_clustalw2(est2genome_fasta):
	clustalw2_string = "clustalw2 -INFILE=" + est2genome_fasta + " -align -type=DNA -outfile=" + curr_align_out \
	+ " -output=phylip-relaxed"
	subprocess(clustalw2_string,shell=True)

def batch_align(refseq_fasta_file,output_dir,velvet_file_dict, number_procs):
	
	align_commands = []
	refseq_fasta = SeqIO.parse(refseq_fasta_file,"fasta")
	pool = mp.Pool(processes=number_procs)
	est2genome_files = [pool.apply( gather_est2genome_hits, args=(refseq,output_dir,velvet_file_dict[refseq.id])) for refseq in refseq_fasta]
	refseq_fasta.close()
	
	results = [pool.apply( run_clustalw2, args=(est2genome_fasta_file)) for est2genome_fasta_file in est2genome_files]




