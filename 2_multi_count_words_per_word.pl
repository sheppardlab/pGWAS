#################################################
# 2nd genotype-phenotype correlation script     #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/local/bin/perl
# give as arguments a folder with contigs and a folder with fasta files containing words
#this script will call the count words script for every word file


## load perl packages
use strict;
use warnings;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Cwd 'abs_path';
use List::MoreUtils qw(firstidx);


sub usage{
	print "---------------------------------------------------------------------------\n";
	
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

sub check_folders{
	my @input_files = @_;
	
	foreach (@input_files){
		if (-d $_){
			#print "Directory'$_' Exists!\n";
		}
		else{
			print "Directory '$_' does not exist!\n";
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

sub print_line{
	for (my $i=0;$i<scalar @_;$i++){
		print $_[$i]."\t";
	}
	print "\n";
}



###################################
###########MAIN####################
###################################

#Total number of arguments must be 2
my $total = $#ARGV + 1;
if ($total != 3) {usage;}

#check that inpt args exist
my @input_files = @ARGV;
check_folders($input_files[0]); # the contigs folder
check_folders($input_files[1]); # the words folder

#create output folder output folder
my $output_folder = $input_files[2];
if (-d $output_folder) {} # output directory directory  exists
else{    # output directory directory  does not exists - create it
	unless(mkdir $output_folder) {
		die "Unable to create $output_folder";
	}
}

#get the directory of the fasta files
my $contigs_directory = abs_path($input_files[0]);
my $words_directory = abs_path($input_files[1]);


#open words folder
opendir (DIR, $words_directory) or die $!;
while (my $gene_words_file = readdir(DIR)) { #for every splited words file file 

	# Use a regular expression to ignore non fasta files
    next if ($gene_words_file =~ m/^\./ || ($gene_words_file !~ m/\.fas$/ && $gene_words_file !~ m/\.fasta$/));

	
	my $cmd = "bjobs | wc -l";
	my  $output = `$cmd`;
	while($output>600){
			$output = `$cmd`;
			sleep(1);
	}	
	
	my @temp = split (/\./,$gene_words_file);
	my $gene_name = $temp[0];
	my $file_path  = $words_directory."/".$gene_words_file;
	
		
	print "calculating $gene_name...\n";	
	
	my $command = "bsub <<< \"perl ../count_words_per_word.pl ".$contigs_directory." ".$file_path." ".$output_folder."\"";
	#print "$command\n";
	system($command);




}
