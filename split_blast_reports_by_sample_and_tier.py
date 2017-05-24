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

	report_handle = open("identify_threshold_report.txt",'w')
	for sample in samples:
		gene_list_lower_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"lower",RefSeqProt,blast_report_dir)
		report_blast_statistics(gene_list_lower_quartile, sample, blast_report_suffix, blast_report_dir, "lower",report_handle)		
		gene_list_median_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"median",RefSeqProt,blast_report_dir)
		
		report_blast_statistics(gene_list_median_quartile, sample, blast_report_suffix, blast_report_dir, "median",report_handle)		

		gene_list_upper_quartile = sample_n_genes_quartile(sample_size,sample,blast_report_suffix,"upper",RefSeqProt,blast_report_dir)

		report_blast_statistics(gene_list_upper_quartile, sample, blast_report_suffix, blast_report_dir, "upper",report_handle)	
	report_handle.close()
	
def report_blast_statistics(gene_hash, sample, blast_report_suffix, blast_report_dir, quartile, report_handle):
	
	for gene in gene_hash.keys():
		blastReport = blast_report_dir + "/" + gene + blast_report_suffix
		report_line = compute_summary_for_gene_and_sample(blastReport, sample, gene) + "\t" + quartile + "\t" + str(gene_hash[gene]) + "\n"
		report_handle.write(report_line)
				
def compute_summary_for_gene_and_sample(blastReport, sample, gene):
	blast = open(blastReport,'r')

	print "sample is:\t" + sample + "\n"
	print "gnee is:\t" + gene + "\n"	

	curr_gene_blob = {'sum_pcident': 0,'max_pcident': 0,'sum_pcpos':0,'max_pcident':0,'avg_counter':0,}
	for line in blast:
		parts = line.split("\t")
		curr_sample = parts[0]
		curr_gene = parts[1]
		#print "curr_gene is:\t" + curr_gene + "\n"
		if(curr_sample == sample and curr_gene == gene):
			print "\tgene and sample matched so updating blob\n"
			update_blob(curr_gene_blob, line)

	avg_pcident = float(curr_gene_blob['sum_pcident']) / float(curr_gene_blob['avg_counter'])
	avg_pcpos = float(curr_gene_blob['sum_pcpos']) / float(curr_gene_blob['avg_counter'])
	report_line = sample + "\t" + gene + "\t" + str(avg_pcident) + "\t" + str(curr_gene_blob['max_pcident'])
	report_line = report_line + "\t" + str(avg_pcpos) + "\t" + str(curr_gene_blob['max_pcpos']) + "\t"

	return report_line

def update_blob(curr_gene_blob, line):

	print "in update, what is blob:\n"
	print curr_gene_blob	
	
	parts = line.split("\t")
	pcident = parts[10]
	pcpos = parts[11]

	if(curr_gene_blob['sum_pcident'] == 0):
		curr_gene_blob['sum_pcident'] = pcident
		curr_gene_blob['sum_pcpos'] = pcpos
		curr_gene_blob['max_pcident'] = pcident
		curr_gene_blob['max_pcpos'] = pcpos
		curr_gene_blob['avg_counter'] = 1
	else:
		curr_gene_blob['sum_pcident'] = curr_gene_blob['sum_pcident'] + pcident
		curr_gene_blob['sum_pcpos'] = curr_gene_blob['sum_pcpos'] + pcpos
		
		if(pcident > curr_gene_blob['max_pcident']):
			curr_gene_blob['max_pcident'] = pcident
		
		if(pcpos > curr_gene_blob['max_pcpos']):
			curr_gene_blob['max_pcpos'] = pcpos
		
		curr_gene_blob['avg_counter'] = curr_gene_blob['avg_counter'] + 1

	
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
	if(len(RefSeq_to_percent_covered_hash.keys()) < sample_size):
		sample_size = len(RefSeq_to_percent_covered_hash.keys())

	print "this many genes:\t" + str(len(tmp_subset_hash.keys())) + "\n"
	print "this many to sample:\t" + str(sample_size) + "\n"
	rand_ints = randint(0,len(tmp_subset_hash.keys()),size=sample_size)
	sampled_hash = {}
	for curr_rand in rand_ints:
		curr_sampled_gene = tmp_subset_hash.keys()[curr_rand]
		curr_sampled_gene_score = tmp_subset_hash[curr_sampled_gene]
		sampled_hash[curr_sampled_gene] = curr_sampled_gene_score
	
	print "sample is:\t" + sample		
	print "sampled_these:\t"
	print sampled_hash

	return sampled_hash	
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
				curr_gene = parts[1]
				if(curr_sample == sample and RefSeq.name == curr_gene):
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


