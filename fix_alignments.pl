#!/usr/bin/perl
#
#Daniel Ence
#06/22/2015

use strict;

my ($RefSeq_file, $alignment_dir) = @ARGV;

my %refseq_hash;
open(RefSeq, $RefSeq_file);
while(<RefSeq>){
	chomp;
	if(/\>([\S]+)\s/){
		$refseq_hash{$1} = 1;		
	}
}

my $header_line="";
foreach my $seqID (keys %refseq_hash){

	my %aln_hash;
	my $header_line="";
	my $sample_num=0;
	my $aln_length=0;
	my $alns_found=0;
	my $curr_aln_pos=0;
	my $ref_seq_pos=0;
	my @aln_pos;	
	my $aln_file = $alignment_dir . "/" . $seqID . ".pptxt";
	open(ALN, $aln_file) or die("Couldn't open this file:\t" . $aln_file);
	while(<ALN>){	
		if(length($header_line) == 0){
			$header_line = $_;
			(my $dum, $sample_num, $aln_length) = split(/\s+/, $header_line);
			
		}								
		else{	
			if($alns_found < $sample_num){
				my $tmp_ID = substr($_,0,10);	
				$aln_pos[$curr_aln_pos]=$tmp_ID;
				my $seq = substr($_,10);
				$seq =~ s/[\s\n]+//g;
				$aln_hash{$tmp_ID} = $aln_hash{$tmp_ID} . $seq;
			
				#Need to remember the position of the refseqi
				# in the alignment
				if($seqID =~ /$tmp_ID/){
					$ref_seq_pos = $curr_aln_pos;	
				}				
		
				$alns_found++;
			}
			else{
				my $seq = $_;
				$seq =~ s/[\s\n]+//g;
				$aln_hash{$aln_pos[$curr_aln_pos]} = 
				$aln_hash{$aln_pos[$curr_aln_pos]} . $seq;
					
			}
			$curr_aln_pos++								
		}
		if($curr_aln_pos == $sample_num){$curr_aln_pos = 0;}
	}

	#After reading the alignment file, need to find the the start/stop codons 
	#in the refSeqalignment
	
	#find the start codon in the refseq
	$aln_hash{$aln_pos[$ref_seq_pos]} =~ /ATG/;
	my $tmp=$aln_pos[$ref_seq_pos];
	print "\n\n\n*****$tmp****\n";
	my $tmp=$aln_hash{$aln_pos[$ref_seq_pos]};
	#print "$tmp\n";
	my $start_pos=$-[0];

	print "###########$seqID##############\n";
	foreach my $sample (keys %aln_hash){
		print "$start_pos\t10\n";
		my $tmp=length($aln_hash{$sample});
		print "$tmp\n";	
		my $tmp= substr($aln_hash{$sample},$start_pos - 10, 10+10);
	#	print "$tmp\n";
		my $tmp = substr($aln_hash{$sample},$start_pos,10);
		print "$sample\t$tmp\n";	
	}	
	my $tmp=1;		
}	
close RefSeq;





