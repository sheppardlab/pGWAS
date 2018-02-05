#################################################
# First genotype-phenotype correlation script   #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/local/bin/perl

#This script takes as an input one or more xmfa files. 
#It creates a folder named xmfa_splited_alignments
#It splits the xmfa files to multiple xfma files 
#Each of these files contains the records that corespond to one gene
#The spliting is done based on the part of the header tha follows after the + symbol
#That part should always be the name of the genes and it must be the same for all the input files
#That way a folder with 1 xmfa file per gene is created

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

sub wrong_header_format{
	print "---------------------------------------------------------------------------\n";
	print "Wrong format in the description of the fasta headers.\n\n";
	print "Fasta headrs should be:\n\t >BIGS_id|gene_alias:startin_pos-endinf_pos + gene_discription\n\n";
	print "And gene discription should be one word and it sould be the same in all the xmfa files.\n";
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

sub print_line{
	for (my $i=0;$i<scalar @_;$i++){
		print $_[$i]."\t";
	}
	print "\n";
}

sub parse_xmfa_header{
	my $header = $_[0];
	my @header_fields = split(/:/, $header);
	my @header_prefix = split(/\|/,$header_fields[0]);
	
	#Leave one of the two uncommented. The first one is printing as a fasta header id:start_pos - end_pos while the second is printing only the id
	#my $parsed_header = $header_prefix[0].":".$header_fields[1];
	my $parsed_header = $header_prefix[0];
	
	return 	$parsed_header;
}

###################################
###########MAIN####################
###################################

#Total number of arguments must be  1 
my $total = $#ARGV + 1;
if ($total != 1) {usage;}

#check if input files exist
my $words_fasta_file = $ARGV[0];
check_files($words_fasta_file);


#create output folder output folder
my $output_folder = "./splited_words/";
if (-d $output_folder) {} # output directory directory  exists
else{    # output directory directory  does not exists - create it
	unless(mkdir $output_folder) {
		die "Unable to create $output_folder";
	}
}


my $current_word_header = "";
my $current_gene = "";
my $previous_gene = "";
my $output_file = "";
my $out;
#open the xmfa file
my $word_file_handler = Bio::SeqIO->new(-file => "$words_fasta_file", -format => 'fasta');

while (my $fasta_record = $word_file_handler->next_seq){

	$current_word_header = $fasta_record->display_id;#get the fasta header
	my @fasta_header = split(/_/, $current_word_header); #split the header
	delete $fasta_header[$#fasta_header]; #delete the last entry	
	$current_gene = join('_',@fasta_header); #join the header which now is the gene name
	
	if($current_gene ne $previous_gene){ #when we have a new gene
		$previous_gene = $current_gene;
		$output_file = ">>$output_folder$current_gene".'.fas';
		$out = Bio::SeqIO->new(-file => $output_file, -format => 'fasta' );
		
	}
	$out->write_seq($fasta_record);
}






