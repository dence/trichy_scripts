#!/home/dence/applications/anaconda/bin/python
#Daniel Ence
#06/27/2015

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import argparse
import re
import subprocess

parser = argparse.ArgumentParser(description="annotate contigs with reference sequences")
parser.add_argument("RefSeq",type=str, help="fasta file of the reference sequences")
parser.add_argument("output",type=str, help="directory to be made and will contain the maker output directories") 
parser.add_argument("directories",nargs='+',type=str, help="list of directories containing velvet contigs")
#take a RefSeq and a bunch of assembly directories
args=parser.parse_args()

print "making directory for output " + args.output
if not os.path.exists(args.output):
	os.mkdir(args.output)

#print args.directories

#figure out which files go with which refseqs
refseq_fasta = SeqIO.parse(args.RefSeq,"fasta") 
filesPerSeqID = {}
samples = []
for refseq in refseq_fasta:
	print "Gathering fasta files of the velvet-assembled contigs"
	new_file = args.output + "/" + refseq.id + ".velvet_contigs.fa"
	new_filehandle = open(new_file, 'w')
	
	filesPerSeqID[refseq.id] = new_file	
	curr_ID_pattern = re.compile(refseq.id)
	#get the velvet contig for reach sample

	for directory in args.directories:
		curr_fasta = directory + "/" + refseq.id + ".classified_reads.fastq_velvet/contigs.fa"
		curr_sample = re.split("/",directory)[1]
		samples.append(curr_sample)
		if(os.path.exists(curr_fasta)):
			curr_fasta_file = SeqIO.parse(curr_fasta,"fasta")
			counter = 0
			for contig in curr_fasta_file:
				curr_id = curr_sample + "_" + str(counter) + "." + refseq.id
				counter = counter + 1
				tmp_record = SeqRecord(contig.seq,id=curr_id,name=curr_id,description="")
				SeqIO.write(tmp_record,new_filehandle,"fasta")
		else:
			print "This file was not found:\t" + curr_fasta + " in this directory:\t" + directory
			continue
	new_filehandle.close()		
	print "Finished gathering fasta files of velvet-assembled contigs"
	#move into the output directory
	subprocess.call("cd " + args.output,shell=True)
	#Run maker for the files that we just made
	subprocess.call("maker -CTL",shell=True)
	#fix the control files for the current refseq
	opts_handle = open("maker_opts.ctl",'r')
	tmp_opts = "tmp_opts"
	tmp_opts_handle = open("tmp_opts",'w')
	
	for line in opts_handle:
		
		if re.match("est=",line):
			tmp_opts_handle.write(line.replace("est= ","est=" + args.RefSeq + " "))
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

	subprocess.call("maker -g " + new_file + " " + tmp_opts + " maker_bopts.ctl maker_exe.ctl",shell=True)

print "Finished running maker on velvet contigs"	
refseq_fasta.close()

refseq_fasta = SeqIO.parse(args.RefSeq,"fasta")
print "Aligning est2genome hits with clustalw"
for refseq in refseq_fasta:
	print "\tAligning for this ref seq:\t"
	curr_datastore = refseq.id + ".velvet_contigs.maker.output/" + refseq.id \
	+ ".velvet_contigs_datastore"
	curr_datastore_log = refseq.id + ".velvet_contigs.maker.output/" + refseq.id \
	+ ".velvet_contigs_master_datastore_index.log" 
	
	est2genome_fasta = args.output + "/" + refseq.id + ".est2genome.fasta"
	est2genome_handle = open(est2genome_fasta,'w')  
	SeqIO.write(refseq,est2genome_handle,"fasta")	

	log_handle = open(curr_datastore_log,'r')
	for line in log_handle:
		seqname = line.split("\t")[0]
		if(re.search(refseq.id,seqname) and re.search("FINISHED",line)):
			print "\n\n" + refseq.id
			print len(refseq.seq)
			seq_dir = line.split("\t")[1]
			#print seq_dir.split("/")
			tmp_refseq = seq_dir.split("/")[3].replace(".","%2E")#hardcoded in this position
			gff_file = refseq.id + ".velvet_contigs.maker.output/" + seq_dir + "/" + tmp_refseq + ".gff"
			gff_handle = open(gff_file,'r')
			for gff_line in gff_handle:
				if(re.search("est2genome",gff_line) and \
				re.search("\texpressed_sequence_match\t",gff_line)):
					curr_start = int(gff_line.split("\t")[3])
					curr_stop = int(gff_line.split("\t")[4])
					curr_strand = gff_line.split("\t")[6]
					should_be_length = curr_stop - curr_start
					
					
					print gff_line
					print should_be_length
	
					tmp_handle = open(filesPerSeqID[refseq.id],'r')	
					tmp_fasta = SeqIO.to_dict(SeqIO.parse(tmp_handle,"fasta")) 
					tmp_handle.close()
					
					if seq_dir.split("/")[3] in tmp_fasta:
						curr_record = tmp_fasta[seq_dir.split("/")[3]]
					else:
						continue
					new_seq = curr_record.seq[curr_start - 1:curr_stop]
					if(curr_strand == "-"):
						new_seq = curr_record.seq[curr_start - 1:curr_stop].reverse_complement()
					print len(new_seq)
					new_record = SeqRecord(new_seq,id=seqname,name=seqname,description="")
					SeqIO.write(new_record,est2genome_handle,"fasta")
			gff_handle.close()
	est2genome_handle.close()	
	
	curr_seqID_align_out = args.output + "/" + refseq.id + ".pptx"
	clustalw2_command = "clustalw2 -INFILE=" + est2genome_fasta + " -align -type=DNA -outfile=" + curr_seqID_align_out \
	+ " -output=phylip-relaxed"
	print "aligning this file:\t" + est2genome_fasta
	print clustalw2_command
	subprocess.call(clustalw2_command,shell=True)	
refseq_fasta.close()
#make a maker directory for curr refseq  
#get the velvet_contig for curr refseq from each assembly directory
#put all the contigs into a fasta file in the maker directory
#make control files with the refseq file as EST evidence and est2genome turned on 
#run maker with the contigs as the genome and the refeseq files as the EST evidence
#get the sequences with their names

