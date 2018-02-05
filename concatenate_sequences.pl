#################################################
# GWAS input file sequence concatenation script #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/local/bin/perl

## load perl packages
use strict;
use warnings;
use Cwd 'abs_path';
use List::MoreUtils qw(firstidx);


sub usage{
	print "---------------------------------------------------------------------------\n\n";
	print "This script automaticaly concatenates a list of dna sequences to one sequence entry.\n";
	print "It takes as an anrgument the file with the sequences (no headers) and the header of the fasta entry.\n";
	print "During execution redirect the output to a new text file.\n\n";
	print "---------------------------------------------------------------------------\n";
	exit;
}

sub check_files{
	my @input_files = @_;
	
	foreach (@input_files){
		if (-f $_){
			#print "File '$_' Exists!\n";
		}
		else{
			print "File '$_' does not exist!\n";
			usage;
		}
	}

}

sub open_file{
	my @input_files = @_;
	
	open(my $fh, "+<:encoding(UTF-8)", $input_files[0]);
	#	|| die "can't open UTF-8 encoded filename: $!";
	return $fh;

}


###################################
###########MAIN####################
###################################

my $total = $#ARGV + 1;
if ($total != 2) {usage;}

check_files($ARGV[0]);
my $fasta_header = $ARGV[1];

#print the fasta header given as an argument
print ">$fasta_header\n";

#open the file with the fasta entries. 
my $file_handler = open_file($ARGV[0]);

my $line = "";
while(<$file_handler>){
$line = $_;
chomp;
$line =~ s/\r\n//g;
$line =~ s/\n//g;
print $line;
}

close $file_handler;

