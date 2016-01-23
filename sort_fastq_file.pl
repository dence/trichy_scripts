#!/usr/bin/perl
#
#Daniel Ence
#January 23, 2016
#

use strict;

my (@fastq_files) = @ARGV;

foreach my $fastq_file (@fastq_files){
	sort_by_readID($fastq_file);
}


###################SUBROUTINES##################
sub sort_by_readID{
	my ($fastq_file) = @_;
	open(FASTQ, $fastq_file);

	my @concat_lines;	
	while(<FASTQ>){
		my $read_ID=$_;
		my $read_BPs=<FASTQ>;
		my $plus=<FASTQ>;
		my $read_qual=<FASTQ>;
		my $concat_line = $read_ID . $read_BPs . $plus . $read_qual;
		#print "concat_line is:\n$concat_line\n";
		push @concat_lines, $concat_line;	
	}
	close FASTQ;

	my @sorted_lines = sort @concat_lines;

	my $sorted_fastq = "sorted." . $fastq_file; 
	open(SORTED_FASTQ, ">$sorted_fastq");	
	foreach my $sorted_line (@sorted_lines){
		#my @curr_lines = split("\n", $sorted_line);
		my ($read_ID, $read_BPs, $plus, $read_qual) = split("\n", $sorted_line);		
		
		print SORTED_FASTQ "$read_ID\n";
		print SORTED_FASTQ "$read_BPs\n";
		print SORTED_FASTQ "$plus\n";
		print SORTED_FASTQ "$read_qual\n";
		
		#foreach my $line (@curr_lines){
		#	print SORTED_FASTQ $line;
		#}
			
	} 

}
