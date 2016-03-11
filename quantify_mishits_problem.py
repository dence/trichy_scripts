#Daniel Ence


from Bio import SeqIO
import sys
import argparse
import os
import errno

parser = argparse.ArgumentParser(description="Report the AA positions covered by at least one HSP for each reference sequence from each sample")
parser.add_argument("--RefSeq",type=str,help="File of Sequences to get results for")
parser.add_argument("--blast_report_suffix", type=str,help="suffix to append to the name of each RefSeq to find the blast report")
args=parser.parse_args()


#open Refseq, make ID->seqlength dictionary
RefSeqFile = SeqIO.parse(args.RefSeq,"fasta")
for RefSeq in RefSeqFile:

	blastReport = RefSeq.name + args.blast_report_suffix
	curr_refseq_length = len(RefSeq.seq)
	AA_position_hash = {}
	print "Refseq\t" + RefSeq.name + "\tis this long:\t" + str(curr_refseq_length) + "\n"
	try:
		#open blastn report
		blastReport = RefSeq.name + args.blast_report_suffix
		currReport = open(blastReport,'r')
		default_hash = {}
		for x in range(0,curr_refseq_length):
			default_hash[x] = 0

		for line in currReport:
			parts = line.split("\t")
			curr_sample = parts[0]	

			print "\tsample is:\t" + curr_sample
			AA_position_hash.setdefault(curr_sample, default_hash)	
			#record which positions were covered by an HSP	
			subject_start_pos = parts[20]
			subject_end_pos = parts[21]
			print "\t\tstarts at\t" + subject_start_pos + "\t" + subject_end_pos + "\n"
			for x in range(int(subject_start_pos),int(subject_end_pos)):
				AA_position_hash[curr_sample][x] = 1
	except EnvironmentError as e: 
	#except:
		print(os.strerror(e.errno))				
		sys.stderr.write("file not found:\t" + blastReport + "\n")
	#compute how much of the AA_positions were covered by an HSP	
	for sample in AA_position_hash:
		curr_sample_coverage_counter = 0
		for x in range(0,curr_refseq_length):
			if(AA_position_hash[sample][x] == 1):
				curr_sample_coverage_counter = curr_sample_coverage_counter + 1
		percent_covered = curr_sample_coverage_counter / curr_refseq_length
		sys.stdout.write(RefSeq.name + "\t" + sample + "\t" + str(percent_covered) + "\n")
