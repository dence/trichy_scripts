#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  To split a large fasta sequence in smaller fasta sequences, with or without sliding window
##########################################################
#depth setting added by Daniel Ence, May 23, 2016
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $version = "1.1";
my $scriptname = "fasta_split.pl";
my $changelog = "
#	- v1.0 = 21 Oct 2015
#	- v1.1 = 26 Oct 2015
			 add -split option
\n";

my $usage = "\nUsage [v$version]: 
    perl $scriptname -in <input.fa> [-len <X>] [-out <output.fa>] [-slide <X>] [-keep <X>] [-uc] [-Nrem] [-v] [-chlog] [h]
	
    PURPOSE:
     Split large fasta sequence(s) from a file in smaller fasta sequences (in a file)
     Several options possible; but if only some sequences are keep, coordinates will be off
     (for example if -keep lc, then chr1_start-end in the output file 
     will correspond to start and end in the \"lowercase genome\")
	  	
    MANDATORY ARGUMENTS:	
     -in <X>  =>   (STRING) input file = fasta file  

    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
     -len     =>   (INT)    length of the split sequences
                            default = 1 kb (1000 nt)
     -out     =>   (STRING) to set output file name
                            default = <in>.split.<len>.fa	
     -slide   =>   (INT)    sliding window of -len nt, insteadt of split the sequences
                            -slide 1 to get all possible sequences (beware of size of files...)   
     -keep    =>   (STRING) chose between uc or lc. relevant only if genome is softmasked
                            this would allow to keep only masked sequences or only non masked sequences
     -uc      =>   (BOOL)   to put sequences in uppercase
                            (can be necessary if lowercases are extracted)
     -Nrem    =>   (BOOL)   to remove Ns. This is done after splitting, 
                            so some sub sequences might be smaller/missing
     -v       =>   (BOOL)   verbose mode, make the script talks to you
     -v       =>   (BOOL)   print version if only option
     -h,-help =>   (BOOL)   print this usage
     -chlog   =>   (BOOL)   print changelog
     -depth   =>   (INT)    desired depth
\n";    
        
################################################################################
# Get arguments/options, check some of them, print details of the run if verbose chosen
################################################################################
my $l = 1000;
my $depth=0;
my ($in,$o,$slide,$keep,$uc,$Nrem,$help,$v,$chlog);
GetOptions ('in=s' => \$in, 'out=s' => \$o, 'len=s' => \$l, 'slide=s' => \$slide, 'keep=s' => \$keep, 'uc' => \$uc, 'Nrem' => \$Nrem, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v, 'depth=s' => \$depth);

#check step for options
die "\n $scriptname v$version\n\n" if ((! $in) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die $usage if (! $in);
die "\n ERROR: $in does not exist?\n\n" if (! -e $in);

unless ($o) {
	$o = $in;
	$o =~ s/(\.fa|\.fasta|\.fas)$//;
	($o = $o.".Nrem") if ($Nrem);
	($o = $o.".uc") if ($uc);
	($slide)?($o = $o.".split.".$l.".".$slide):($o = $o.".split.".$l);
	($o = $o.".".$keep) if ($keep);
	$o = $o.".fa";
}
#verbose
if ($v) {
	print STDERR "\n--------------------------------------------------\n";
	print STDERR " --- $scriptname started (v$version)\n";
	print STDERR "     Sequence(s) from: $in\n"; 
	print STDERR "     will be split in sequences of: $l nt\n";
	print STDERR "     with a sliding window of: $slide nt\n" if ($slide);	
	print STDERR "     and printed in: $o\n";
	print STDERR "     Also:\n";	
	print STDERR "      - N will be removed\n" if ($Nrem);	
	print STDERR "      - only uppercase sequences will be kept\n" if (($keep) && $keep eq "uc");	
	print STDERR "      - only lowercase sequences will be kept\n" if (($keep) && $keep eq "lc");	
	print STDERR "      - sequences will be rewritten in uppercase letters\n" if ($uc);	
	print STDERR " --- Splitting in progress...\n";
}
$keep = "na" unless ($keep);
$slide = "na" unless ($slide);
($uc)?($uc="y"):($uc="n");
($Nrem)?($Nrem="y"):($Nrem="n");

my $total = fasta_split($in,$l,$keep,$slide,$uc,$Nrem,$v,$depth);

if ($v) {
	print STDERR "     ...done\n";
	print STDERR "     total amount of sequence written = $total nt\n";
	print STDERR "--------------------------------------------------\n\n";
}	
exit;

##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# fasta_split
# fasta_split($in,$l,$keep,$uc,$Nrem,$v);
#----------------------------------------------------------------------------
sub fasta_split {
	($in,$l,$keep,$slide,$uc,$Nrem,$v,$depth) = @_;
	my $total = 0;
	
	#set up the depth vs read lengt#set up depth vs read length
	my $percent_seqed = 0;
	my $times_seqed = 0;

	#print "Depth is:\t" . $depth . "\n";

	if($depth > $l){
		#if depth > read length, means doing some seqs twice
		#
	        #for example if want 48X with 32bp reads
	        #need to double on 0.5 half of the reads
	        #need to handle the remainder too
	        $times_seqed = int($depth / $l);
		$percent_seqed = $depth / $l - int($depth / $l);
	}
	else{   
	        #if depth < read length 
	        $percent_seqed = $depth / $l;
	}

	#print "percent_seqed is:\t" . $percent_seqed . "\n";
	#print "times_seqed is:\t" . $times_seqed . "\n";
	#print "length is:\t" . $l . "\n";

	my $fa = Bio::SeqIO->new(-file => $in, -format => "fasta") or confess "     \nERROR (main): Failed to create Bio::SeqIO object from $in $!\n";
	open(my $o_fh, ">", $o) or confess "     \nERROR (main): could not open to write $o $!\n";
	SEQ: while( my $seq = $fa->next_seq() ) {
		my $id = $seq->display_id;
		my $desc = $seq->desc;	
		my $fas = $seq->seq;
		$fas =~ s/N//g if ($Nrem eq "y"); #remove Ns if relevant	
		$fas =~ s/[A-Z]//g if ($keep eq "lc"); #remove uc if lc should be kept
		$fas =~ s/[a-z]//g if ($keep eq "uc"); #remove lc if uc should be kept		
		$fas = uc($fas) if ($uc eq "y"); #re write in uc if relevant
		my $fal = length($fas);
		$total+=$fal;
   		my $seqobj = Bio::Seq->new( -display_id => $id, -seq => $fas);
   		my $inc;
   		($slide eq "na")?($inc = $l):($inc = $slide);
		SPLIT: foreach (my $st = 1; $st < $fal; $st+=$inc) {
			my $en = $st+$l-1;
			$en = $fal if ($en > $fal); #to make sure; smaller sequence for the last one

			my $split = $seqobj->subseq($st => $en); #subseq
			
			if($depth == 0){
				my $newid = $id."_".$st."-".$en." ".$desc;
				print $o_fh ">".$newid."\n".$split."\n";
			}
			else{

				if($times_seqed > 0){			

					for(my $i=0; $i<=$times_seqed; $i++){
						#my $newid = $id."_".$st."-".$en." ".$desc;
						my $newid = $id."_".$st."-".$en.".".$i." ".$desc;
						print $o_fh ">".$newid."\n".$split."\n";
					}
				}

				if($percent_seqed > 0){
					
					if(rand(1) < $percent_seqed){
						my $newid = $id."_".$st."-".$en.$percent_seqed." ".$desc;
						print $o_fh ">".$newid."\n".$split."\n";
					}			
				}
			}
			next SEQ if (($slide ne "na") && ($en > $fal)); #no need to slide more in that case
		}
	}
	close ($o_fh);
	return ($total);
}
