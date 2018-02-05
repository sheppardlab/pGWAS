#! /usr/local/bin/perl

## load perl packages
use strict;
use warnings;
use Cwd 'abs_path';
use List::MoreUtils qw(firstidx);


sub usage{
	print "---------------------------------------------------------------------------\n\n";
	print "This script automaticaly parces the supergenome input files and produces a tab seperated file with a list of all genes (with their description and sequence).\n";
	print "It takes as an anrgument only the folder with the files.\n";
	print "During execution redirect the output to a new text file.\n\n";
	print "---------------------------------------------------------------------------\n";
	exit;
}

sub open_file{
	my @input_files = @_;
	
	open(my $fh, "+<:encoding(UTF-8)", $input_files[0]);
	#	|| die "can't open UTF-8 encoded filename: $!";
	return $fh;

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

###################################
###########MAIN####################
###################################

my $total = $#ARGV + 1;
if ($total != 1) {usage;}

check_folders($ARGV[0]);
my $directory = abs_path($ARGV[0]); 

opendir (DIR, $directory) or die $!;


my $current_file = "";
my $current_file_handler = "";
my $line = "";
my $first_line_flag;
my $gene_name;
my $gene_descreption;
my $gene_sequence;
while (my $file = readdir(DIR)) { #for every file in the directory
	
	#skip if the file is not fas or txt
	next if ($file =~ m/^\./ || ($file !~ m/\.fas$/ && $file !~ m/\.txt$/));
	
	
	$current_file = $directory ."/". $file;
	$current_file_handler = open_file($current_file);
	
	$line = "";
	$first_line_flag = 0;
	while(<$current_file_handler>){ #for every line of the file
		$line =$_;
		chomp;
		
		#skip the first line
		if($first_line_flag == 0){
			$first_line_flag = 1; next;
		}
		
		#initialise/clear the values before every line
		$gene_name = "";
		$gene_descreption = "";
		$gene_sequence = "";
		
		my @temp = split('\t', $line);
		$gene_name = $temp[0];
		$gene_descreption = $temp[3];
		$gene_sequence = uc $temp[12];
		
		#print the entry 
		my $line = "$gene_name\t$gene_descreption\t$gene_sequence"; 
		$line =~ s/upodiastoli/_/g;
		print $line; 
		
	
	
	}
	
	close $current_file_handler;
	

}