#################################################
# 3rd genotype-phenotype correlation script     #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/local/bin/perl
#takes as an input a contigs folder, a fasta file with the associated per gene words (teh fasta file should have teh name of the gene)
#and an output folder 
#The output is one file per word with the presence or absence of the specifi word in the strains in the contigs folder


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
check_folders($input_files[0]);
check_files($input_files[1]);

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
my $words_file = $input_files[1];

my $current_fasta_file = ""; #here the name of the fasta file will be saved
my $current_abs_path_fasta_file = ""; #here the absolut path of the fasta file will be saved
my $current_strain_id = "";
my $concatenated_fasta = "";
my $word_counter = 0;



#open the words file 
my $word_seq_io = Bio::SeqIO->new(-file => "$words_file", -format => 'fasta');
while (my $seq = $word_seq_io->next_seq){ #for every_word
	
	my $word = $seq->seq; #print $word;
	my $reverce_complement = $seq->revcom; 
	my $reverce_complement_word = $reverce_complement->seq; 
	my $word_header = $seq->display_id; #print "$word_header\n";
	my $output_file = $output_folder."\/".$word_header.".txt"; #print "$output_file\n";
	
	#create_output_file
	open(my $fh, '>', $output_file) or die "Could not open file '$output_file' $!";
	
	#open contig folder
	opendir (DIR, $contigs_directory) or die $!;
	while (my $file = readdir(DIR)) { #for every contig file 
	
		# Use a regular expression to ignore non fasta files
		next if ($file =~ m/^\./ || ($file !~ m/\.fas$/ && $file !~ m/\.fasta$/));
		
		#clear the variables
		$current_fasta_file = "";
		$current_abs_path_fasta_file = "";
		$current_strain_id = "";
		$concatenated_fasta = "";
		$word_counter = 0;
		
		$current_fasta_file = $file;
		$current_abs_path_fasta_file = $contigs_directory ."/". $current_fasta_file;
		my @temp = split('\.', $file);
		$current_strain_id = $temp[0];

		#for every sequence in the fasta file concatenate the contigs
		my $seq_io = Bio::SeqIO->new(-file => "$current_abs_path_fasta_file", -format => 'fasta');
		while (my $seq = $seq_io->next_seq){
			my $qurrent_entry = $seq->seq;
			$concatenated_fasta .= $qurrent_entry; #the concatenated sequence of each isolate will be stored here
		}
		
		if ($concatenated_fasta =~ /$word/){
			$word_counter ++; 
			#print "found in $file\n$concatenated_fasta\n$word\n"; exit;
		}
		elsif ($concatenated_fasta =~ /$reverce_complement_word/){
			$word_counter ++;
			#print "found RC in $file\n$concatenated_fasta\n$word\n"; exit;
		}
	
		print $fh "$current_strain_id\t$word_counter\n";
		
	
	
	}
	close $fh;
		
}






