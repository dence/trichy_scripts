#!/usr/bin/perl

#Daniel Ence November 2, 2015


use strict;

####################### MAIN###################################

my %file_prefix_hash= (

"Pentatrich.hominis"=>"../sven_gould_trich_species_data/cleaned_fastqs/Pentatrich.homi.PhGII.gould.qc.fastq",
"Tetratrich.gall"=>"../sven_gould_trich_species_data/cleaned_fastqs/Tetratrich.gall.M3.gould.qc.fastq",
"Trich.gall"=>"../sven_gould_trich_species_data/cleaned_fastqs/Trich.gall.GCB.gould.qc.fastq",
"Trich.tenax"=>"../sven_gould_trich_species_data/cleaned_fastqs/Trich.tenax.HS-4.gould.qc.fastq",
"Trichomitus.batrachorum"=>"../sven_gould_trich_species_data/cleaned_fastqs/Trichomitus.batr.BUB.gould.qc.fastq",

"Tvag.G3"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/G3_Tv.prinseq.qc_trimmed.fastq",
"Tvag.30235"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/30235_Tv.prinseq.qc_trimmed.fastq",
"Tvag.B7RC2"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/B7RC2_Tv.prinseq.qc_trimmed.fastq",
"Tvag.C1"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/C1_Tv.prinseq.qc_trimmed.fastq",
"Tvag.JRSTV41"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/JRSTV41_Tv.prinseq.qc_trimmed.fastq",
"Tvag.T1"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/T1_Tv.prinseq.qc_trimmed.fastq",

"Tvag.30001"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.vag.30001.prinseq.qc_trimmed.fastq",
"Tvag.30238"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.vag.30238.prinseq.qc_trimmed.fastq",
"Tvag.50143"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.vag.50143.prinseq.qc_trimmed.fastq",
"Tvag.t016"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.vag.t016.prinseq.qc_trimmed.fastq",
"T.foetus"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.foetus.prinseq.qc_trimmed.fastq",
"T.hominis"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.hominis.prinseq.qc_trimmed.fastq",
"T.tenax"=>"../prinseq_fqtrimmed_fastqs/second_batch_samples/T.tenax.prinseq.qc_trimmed.fastq"

);

#my %file_prefix_hash = (
#
#"Tvag.G3"=>"../prinseq_fqtrimmed_fastqs/first_batch_samples/reads_for_this_repark/G3_Tv.prinseq.qc_trimmed.fastq"
#
#);


my @kmer_lengths = (30,27,24,21);



#foreach sample
foreach  my $sample (keys %file_prefix_hash){
	print "key is:\t" . $sample . "\n";
	my $val = %file_prefix_hash->{$sample};
	print "value is:\t" . $val . "\n"; 
#foreach kmer length
	foreach my $kmer_length (@kmer_lengths){

		my $filepath = %file_prefix_hash->{$sample};
		#make the config file for the sample and kmer length
		my $curr_config_file = make_config_file($sample, $filepath,$kmer_length);	
		#run taxnomer for the sample and kmer length
		my $kmer_dir = "single_CC_prot_k" . $kmer_length;
		print "cd to $kmer_dir\n"; 
		chdir($kmer_dir);
		system("classify_reads " . $curr_config_file);
		chdir("..");	
	}	
}

######################### SUBROUTINES #######################

sub make_config_file{
	
	my ($sample, $file_path, $kmer_length) = @_;

	open(TEMPLATE,"taxonomer.config.protein.template");
	my $kmer_folder="single_CC_prot_k" . $kmer_length;
	my $new_config_file ="taxonomer.config." . $sample . ".k" . $kmer_length;
	my $new_config_path = $kmer_folder . "/" . "taxonomer.config." . $sample . ".k" . $kmer_length;
	open(NEW_CONFIG, ">$new_config_path") or die "Couldn't open this file for writing:\t$new_config_path\n";
	print "opened this file for writing:\t" . $new_config_file . "\n";
	while(<TEMPLATE>){
		#print "sending this to substitute_arg:\t" . $_;
		#print "$sample\n";
		#print "$file_path\n";
		#print "$kmer_length\n";
		my $new_line = substitute_arg($_,$sample,$file_path, $kmer_length);
		#print "NEW LINE IS:\t" . $new_line;
		print NEW_CONFIG $new_line;	
	}
	close TEMPLATE;
	close NEW_CONFIG;
	print "finished writing to this file:\t" . $new_config_file . "\n";
	return $new_config_file;
}

sub substitute_arg{
	my ($line, $sample, $fastq_path, $kmer_length) = @_;

	my $new_line = $line;
	if($new_line =~ /dbXX/){
		my $currDB="Single_CC_Library.prot.k" . $kmer_length;
		$new_line =~ s/dbXX/$currDB/;
	}
	elsif($new_line =~ /stiXX/){
		my $currSTI="../Single_CC_Library.prot.sti";
		$new_line =~ s/stiXX/$currSTI/;
	}
	elsif($new_line =~ /triXX/){
		my $currTRI = "../Single_CC_Library.prot.tri";
		$new_line =~ s/triXX/$currTRI/;
	}
	elsif($new_line =~ /kcXX/){
		my $currKC = "Single_CC_Library.prot_taxonomer.k" . $kmer_length . ".kc";
		$new_line =~ s/kcXX/$currKC/;
	}
	elsif($new_line =~ /inXX/){
		my $currIN = "../Single_CC_Library.prot_taxonomer.fa";
		$new_line =~ s/inXX/$currIN/;
	}
	elsif($new_line =~ /fastqXX/){
		my $currFastq = "../" . $fastq_path;
		$new_line =~ s/fastqXX/$currFastq/;
	}
	elsif($new_line =~ /outXX/){
		my $currOUT = $sample . ".k" . $kmer_length . ".taxonomer.out";
		$new_line =~ s/outXX/$currOUT/;
	}
	
	return $new_line;	
}	

