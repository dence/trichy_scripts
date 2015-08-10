from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import re
import subprocess

import multiprocessing as mp

import pdb


def write_contigs(to_filehandle,from_filehandle,curr_sample,curr_refseq_id):

	counter = 0
	for contig in from_filehandle:
		curr_id = curr_sample + "_" + str(counter) + "." + curr_refseq_id
		counter = counter + 1
		tmp_record = SeqRecord(contig.seq,id=curr_id,name=curr_id,description="")
		SeqIO.write(tmp_record,to_filehandle,"fasta")

def gather_refseq_contigs(refseq_obj,output_dir,directories):
	
	new_fasta_file = output_dir + "/" + refseq_obj.id + ".velvet_contigs.fa"
	new_filehandle = open(new_fasta_file, 'a')

	print "Gathering to this fasta file:\t" + new_fasta_file

	wrote_Bool = 0
	for curr_dir in directories:
		curr_velvet_file = curr_dir + "/" + refseq_obj.id + ".classified_reads.fastq_velvet/contigs.fa"
		curr_sample = re.split("/",curr_dir)[len(re.split("/",curr_dir)) - 1]
		if(os.path.exists(curr_velvet_file)):
			
			#can't use the dictionary method because 
			#we sometimes get multiple contigs
			#curr_fasta_dict = SeqIO.to_dict(SeqIO.parse(curr_velvet_file,"fasta"))
			curr_fasta_file = SeqIO.parse(curr_velvet_file,"fasta")
			write_contigs(new_filehandle,curr_fasta_file,curr_sample,refseq_obj.id)
			wrote_Bool = 1
		else:
			print "This file was not found:\t" + curr_velvet_file + " in this directory:\t" + curr_dir
	new_filehandle.close()

	if(wrote_Bool == 1):
		return refseq_obj.id + ".velvet_contigs.fa" 

def gather_and_annotate_refseq_contigs(arg_duples):

	refseq_obj = arg_duples[0]
	output_dir = arg_duples[1]
	directories = arg_duples[2]

	print "Starting on this refseq:\t" + refseq_obj.id

	#This gathers the contigs that were assembled by velvet for that refseq
	velvet_fasta_file = gather_refseq_contigs(refseq_obj, output_dir, directories)
	

	if(velvet_fasta_file):
		#This annotates the contigs that we just gathered	
		refseq_file = arg_duples[3]
		annotate_velvet_contigs(velvet_fasta_file)
		gather_and_align_est2genome_hits((refseq_obj, output_dir, velvet_fasta_file))	
	print "Finished with this refseq:\t" + refseq_obj.id

def gather_annotate_and_align_batch(RefSeq_file,output_dir,directories,number_procs):
	
	#figuring out which files go with which refseqs	
	refseq_fasta = SeqIO.parse(RefSeq_file,"fasta")

	command_duples = []
	for refseq in refseq_fasta:
		command_duples.append((refseq, output_dir, directories,RefSeq_file))
	refseq_fasta.close()
	os.chdir(output_dir)

	print "CWD now is:\t" + os.getcwd()

	prepare_opts_file(RefSeq_file)

	pool = mp.Pool(processes=number_procs)
	pool.map(gather_and_annotate_refseq_contigs, command_duples)	
	print "finished annotating velvet_contigs\n"
	
def prepare_opts_file(RefSeq_file):

        subprocess.call("maker -CTL",shell=True)

        opts_handle = open("maker_opts.ctl",'r')
        tmp_opts = "tmp_opts"
        tmp_opts_handle = open(tmp_opts,'w')

        print "preparing the tmp_opts file\n"
        for line in opts_handle:
                if(re.match("est=",line)):
                        tmp_opts_handle.write(line.replace("est= ","est=" + "../" + RefSeq_file + " "))
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

def annotate_velvet_contigs(velvet_fasta_file):	
		
	print "Running maker for:\t" + velvet_fasta_file
	maker_command_string = "maker -g " + velvet_fasta_file + " tmp_opts maker_bopts.ctl maker_exe.ctl"
	subprocess.call(maker_command_string,shell=True)
	
def gather_est2genome_seqs(refseq_obj, est2genome_fasta, log_line, velvet_file):
		
	seq_dir = log_line.split("\t")[1]
	tmp_refseq = seq_dir.split("/")[len(seq_dir.split("/")) - 2].replace(".","%2E")#hardcoded in this position
	gff_file = refseq_obj.id + ".velvet_contigs.maker.output/" + seq_dir + "/" + tmp_refseq + ".gff"
	gff_handle = open(gff_file,'r')
	for gff_line in gff_handle:
		if(re.search("est2genome",gff_line) and \
		re.search("\texpressed_sequence_match\t",gff_line)):
			curr_start = int(gff_line.split("\t")[3])
			curr_stop = int(gff_line.split("\t")[4])
			curr_strand = gff_line.split("\t")[6]
			
			tmp_handle = open(velvet_file,'r')
			tmp_fasta = SeqIO.to_dict(SeqIO.parse(tmp_handle,"fasta"))
			tmp_handle.close()
	
			if seq_dir.split("/")[len(seq_dir.split("/")) - 2] in tmp_fasta:
				curr_record = tmp_fasta[seq_dir.split("/")[len(seq_dir.split("/")) - 2]]
				seqname = curr_record.id
			else:
				continue
			new_seq = curr_record.seq[curr_start - 1:curr_stop]
			if(curr_strand == "-"):
				new_seq = curr_record.seq[curr_start - 1:curr_stop].reverse_complement()
			new_record = SeqRecord(new_seq,id=seqname,name=seqname,description="")
			
			SeqIO.write(new_record, open(est2genome_fasta,'a'), "fasta")
	
def gather_est2genome_hits(refseq_obj, output_dir, velvet_file):
	
	est2genome_fasta = refseq_obj.id + ".est2genome.fasta"
	est2genome_handle = open(est2genome_fasta,'a')

	curr_datastore = output_dir + "/" + refseq_obj.id + ".velvet_contigs.maker.output/" + refseq_obj.id \
	+ ".velvet_contigs_datastore"
	curr_datastore_log = output_dir + "/" + refseq_obj.id + ".velvet_contigs.maker.output/" + refseq_obj.id \
	+ ".velvet_contigs_master_datastore_index.log"
	est2genome_handle.close()
	
	if(os.path.exists(curr_datastore_log)):
		log_handle = open(curr_datastore_log,'r')
		for line in log_handle:
			seqname = line.split("\t")[0]
			if(re.search(refseq_obj.id,seqname) and re.search("FINISHED",line)):
				gather_est2genome_seqs(refseq_obj, est2genome_fasta,line, velvet_file)
		log_handle.close()	
		return est2genome_fasta
	else:
		print "couldn't find this datastore log:\t" + curr_datastore_log

def run_clustalw2(est2genome_fasta, refseq_obj):
	print "Running clustalw2 for:\t" + est2genome_fasta
	curr_align_out = refseq_obj.id + ".est2genome.pptx"
	clustalw2_string = "clustalw2 -INFILE=" + est2genome_fasta + " -align -type=DNA -outfile=" + curr_align_out \
	+ " -output=phylip-relaxed"
	print clustalw2_string
	subprocess.call(clustalw2_string,shell=True)

def gather_and_align_est2genome_hits(arg_duple):

	refseq_obj = arg_duple[0]
	output_dir = arg_duple[1]
	velvet_file = arg_duple[2]
	est2genome_fasta = gather_est2genome_hits(refseq_obj, output_dir, velvet_file)

	
	if(est2genome_fasta):
		run_clustalw2(est2genome_fasta, refseq_obj)
	else: 
		print "no est2genome hits for:\t" + refseq_obj.id

	
