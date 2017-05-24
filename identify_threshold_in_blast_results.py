#Daniel Ence
#March 23, 2016

import sys
import os.path
import re
import argparse
import subprocess
from scipy.stats.mstats import mquantiles
from numpy.random import randint
from Bio import SeqIO

def main(RefSeqProt,RefSeqTrans,blast_report_suffix,list_of_samples,sample_size,blast_report_dir,fastq_filter_dir):

	#get the list of samples out of the list_of_samples_string
	samples = list_of_samples.split(",")

	for sample in samples:
		gene_list_lower_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"lower",RefSeqProt,blast_report_dir)
		align_fastqs_for_gene(sample,gene_list_lower_quartile,fastq_filter_dir,RefSeqTrans,sample + "_" + "lower")
		
		gene_list_median_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"median",RefSeqProt,blast_report_dir)
		align_fastqs_for_gene(sample,gene_list_median_quartile,fastq_filter_dir,RefSeqTrans,sample + "_" + "median")
		
		gene_list_upper_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"upper",RefSeqProt,blast_report_dir)
		align_fastqs_for_gene(sample,gene_list_upper_quartile,fastq_filter_dir,RefSeqTrans,sample + "_" + "upper")
		

def align_fastqs_for_gene(sample,gene_list,fastq_filter_dir,RefSeqFilename,working_dir_name):
	#for each gene in each list
		#get the fastq for that gene for this sample
		#align them to the refseq with a very sensitive bowtie2 alignment
		#record how many of the reads were align-able	
		#build a consensus sequence with pileup
		#blast that consensus against the refseqs
		#report how well the 
	#make a directory to hold the alignments
	working_dir = working_dir_name + "_working_dir"
	if(not os.path.exists(working_dir)):
		os.mkdir(working_dir) 
	for gene in gene_list:
		fastq = fastq_filter_dir + "/" + sample + "/" + gene + ".classified_reads.fastq"
		#print "this is the fastq:\t" + fastq
		if(os.path.exists(fastq)):
		
			bwa_or_bowtie2 = "bwa"
			
			#print "and we found it!\n"
			#make a reference fasta for the alignment
			sequence_dictionary = SeqIO.to_dict(SeqIO.parse(open(RefSeqFilename,'r'),"fasta"))
			tmp_list = [sequence_dictionary[gene]]
			tmp_fasta = working_dir + "/" + gene + ".tmp.fa"
			tmp_fasta_file = open(tmp_fasta,'w')	
			SeqIO.write(tmp_list,tmp_fasta_file,"fasta")
			tmp_fasta_file.close()
			
			align_fastqs(fastq,tmp_fasta, sample, bwa_or_bowtie2)

def align_fastqs(fastq,tmp_fasta_name, sample, bwa_or_bowtie2):

	sam_file = "" 
	if(bwa_or_bowtie2 == "bowtie2"):	
		build_string = "bowtie2-build " + tmp_fasta_name + " " + tmp_fasta_name
		subprocess.call(build_string,shell=True)
		sam_file = tmp_fasta_name + "_" + sample + ".sam"
		align_out = tmp_fasta_name + "_bowtie2.out"
		bowtie2_string = "bowtie2 -x " + tmp_fasta_name + " --local -D 50 -N 1 -R 3 -L 18 -i S,1,0.50 -U " + fastq  + " -S " + sam_file + " &> " + align_out 
		print "this is the bowtie2 command:\t" + bowtie2_string
		subprocess.call(bowtie2_string,shell=True)
		#align the fastq to the tmp_fasta with bowtie2?			
		#build a consensus with pileup
		#need to faidx the fasta file
	elif(bwa_or_bowtie2 == "bwa"):
		build_string = "bwa index " + tmp_fasta_name
		print "build_string is:\t" + build_string
		subprocess.call(build_string,shell=True)
		sam_file = tmp_fasta_name + "_" + sample + ".bwa_mem.sam"
		align_out = tmp_fasta_name + "_bwa_mem.out"
		bwa_string = "bwa mem " + tmp_fasta_name + " " + fastq + " > " + sam_file + " 2> " + align_out
		print "this is the bwa command:\t" + bwa_string
		subprocess.call(bwa_string,shell=True)
	
	 
	faidx_command = "samtools faidx " + tmp_fasta_name
	print faidx_command
	subprocess.call(faidx_command,shell=True)
	bam_file = sam_file + ".bam"
	view_string = "samtools view -bSo " + bam_file + " " + sam_file
	subprocess.call(view_string, shell=True) 
	sorted = sam_file + ".sorted"
	sort_command = "samtools sort " + bam_file + " " + sorted
	print sort_command
	subprocess.call(sort_command,shell=True) 
	consensus_file = tmp_fasta_name + "_" + sample + "_consensus.out"
	consensus_error = tmp_fasta_name + "_" + sample + "_consensus.error"
	
	consensus_command = "samtools mpileup -uf " + tmp_fasta_name + " " + sam_file + " > " + consensus_file + " 2> " + consensus_error 
	#consensus_command = "samtools mpileup -uf " + tmp_fasta_name + " " + sam_file + " | /home/dence/applications/bcftools/bcftools view -cg - | /home/dence/applications/bcftools/vcfutils.pl vcf2fq > " + consensus_file + " 2> " + consensus_error
	print consensus_command + "\n"
	subprocess.call(consensus_command,shell=True) 
						
	
def sample_n_genes_quartile(sample_size, sample, blast_report_suffix,quartile, RefSeq, blast_report_dir):
	#print "sampling this many:\t" + str(sample_size) + " genes for this sample\t" + sample + "\n"

	RefSeq_to_percent_covered_hash = make_refseq_to_percent_covered_hash(RefSeq,blast_report_suffix,blast_report_dir,sample)
	percent_covered_keys = RefSeq_to_percent_covered_hash.keys()
	percent_covered_scores = RefSeq_to_percent_covered_hash.values()
	
	quantiles = mquantiles(percent_covered_scores)
	
	range = ()
	if(quartile == "lower"):
		range = (0,quantiles[0])			
	elif(quartile == "median"):
		range = (quantiles[0],quantiles[1])	
	elif (quartile == "upper"):
		range = (quantiles[1], max(percent_covered_scores))

	print "quartile is:\t" + quartile + "\n"
	print "range is:\n"
	print range 
	tmp_subset_hash = {}
	subset_size = 0
	for gene in RefSeq_to_percent_covered_hash.keys():
		curr_score = RefSeq_to_percent_covered_hash[gene]
		if(curr_score > range[0] and curr_score < range[1]):
			subset_size = subset_size + 1		
			tmp_subset_hash[gene] = RefSeq_to_percent_covered_hash[gene]
			#print "subset_size:\t" + str(subset_size) + "\n"
		
	print "this many genes in range:\n"
	print len(tmp_subset_hash.keys())

	sampled_counter = 0
	sampled_list = {}
	if(len(RefSeq_to_percent_covered_hash.keys()) < sample_size):
		sample_size = len(RefSeq_to_percent_covered_hash.keys())

	print "this many genes:\t" + str(len(tmp_subset_hash.keys())) + "\n"
	print "this many to sample:\t" + str(sample_size) + "\n"
	rand_ints = randint(0,len(tmp_subset_hash.keys()),size=sample_size)
	sampled_list = []
	for curr_rand in rand_ints:
		sampled_list.append(tmp_subset_hash.keys()[curr_rand])
		
	print "sampled_these:\t"
	print sampled_list

	return sampled_list	
	#find the quartile out of the blast reports USE SCIPY!!!
		
	#grab "sample_size" number of genes out of that quartile
	#return a list of those genes	

def make_refseq_to_percent_covered_hash(RefSeqFilename,blast_report_suffix, blast_report_dir,sample):

	RefSeqFile = SeqIO.parse(RefSeqFilename,"fasta")
	percent_covered_hash = {}	
	for RefSeq in RefSeqFile:
		try:
			curr_refseq_length = len(RefSeq.seq)
			blastReport = blast_report_dir + "/" + RefSeq.name + blast_report_suffix
			currReport = open(blastReport,'r')	
			default_hash = {}
			for x in range(0,curr_refseq_length):
				default_hash[x] = 0
			
			for line in currReport:
				parts = line.split("\t")
				curr_sample = parts[0]
				if(curr_sample == sample):
					subject_start_pos = parts[20]
					subject_end_pos = parts[21]

					for x in range(int(subject_start_pos),int(subject_end_pos)):
						default_hash[x] = 1
			curr_sample_coverage_counter = 0		
			for x in range(0,curr_refseq_length):
				if(default_hash[x] == 1):
					curr_sample_coverage_counter = curr_sample_coverage_counter + 1		
			percent_covered = float(curr_sample_coverage_counter) / float(curr_refseq_length)
			if(percent_covered > 0):
				#print "this gene:\t" + RefSeq.name + "\thad this much covered\t" + str(percent_covered) + "\n"
				percent_covered_hash[RefSeq.name] = percent_covered	
		except EnvironmentError as e:
			dumb = 0	
	return percent_covered_hash		

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="select a given number of genes from the lower, median and upper quartiles of blast results, align the reads classified to those genes to the refseqs, build consensus and blast the consensus seqs to the refseqs. Answer questions about how many reads aligned, how well did the consensus seq match the refseq ")
	#parameters
	#blast reports
	#list of samples to analyze
	#number of genes to select
	#reference sequences
	#location of directory with the filtered fastqs for the samples
	parser.add_argument("--RefSeqProt",type=str,help="File of Reference Amino Acid Sequences used for classification")
	parser.add_argument("--RefSeqTrans",type=str,help="File of Reference Transcript Sequences used for classification")
	parser.add_argument("--blast_report_dir",type=str,help="Directory that holds the blast reprots")
	parser.add_argument("--blast_report_suffix",type=str,help="suffix to append to the name of each RefSeq to find the balst report")
	parser.add_argument("--list_of_samples",type=str,help="comma separated list of samples to do this analysis for")
	parser.add_argument("--sample_size",type=int,help="integer to give number of genes to be sampled")
	parser.add_argument("--fastq_filter_dir",type=str,help="directory with the fastq files")
	args=parser.parse_args()

	main(args.RefSeqProt,args.RefSeqTrans,args.blast_report_suffix,args.list_of_samples,args.sample_size,args.blast_report_dir,args.fastq_filter_dir)


