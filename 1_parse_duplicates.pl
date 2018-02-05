#! /usr/local/bin/perl

## load perl packages
use strict;
use warnings;
use Cwd 'abs_path';
use List::MoreUtils qw(firstidx);

#Define those variables before executing the script. 
my $folder_prefix = "supergenome_input/pangenome_union.txt_";
my $reference_file_name = "supergenome_input/rp62a_old.txt_";

#my $reference_genome_gene_name_pattern = "(b\d\d\d\d)_";

sub usage{
	print "---------------------------------------------------------------------------\n\n";
	print "This script automaticaly parces the duplicates.log output file of the supergenome script.\n";
	print "It takes as an onrument only the duplicates file name.\n";
	print "Before execution open the source code and change the patterns in the begenning of the document.\n";
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


check_files($ARGV[0]);

my $duplicates_file = open_file($ARGV[0]);
my $line = "";
my $first_line_flag = 0;

while(<$duplicates_file>){ #for every line
	$line =$_;
	chomp;
 
	#skip the first line
	if($first_line_flag == 0){
		$first_line_flag = 1; next;
	}

	#remove the folder prefix
	$line =~ s/$folder_prefix//g;
	
	#remove the reference file name
	$line =~ s/$reference_file_name//g;
	

	#remove the file name from the beginning of each line
	$line =~ s/00[0-9]+_[0-9]+\.txt_//g;
	
	
	
	#add a tab_after every reference_gene_name
	
	$line =~ s/(SERP[0-9]+)_/$1\t/;
	$line =~ s/(SE[0-9]+)_/$1\t/;
	$line =~ s/(SEA[0-9]+)_/$1\t/;
	$line =~ s/(SEupodiastolip[0-9]+)_/$1\t/;
	#$line =~ s/(ECs[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(pOSAK1upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(EcHSupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(APECO1upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(APECO78upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(ETECupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(N1917upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(EC55989upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(ECS88upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(pECS88upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(ECUMNupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(p1ECUMNupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(p2ECUMNupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(SBOupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(SDYupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(p53638-226upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(p53638-75upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(O127upodiastolipMAR2upodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(UMNK88upodiastolipEntupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(pAPECO1CoBMupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(pAPECO1Rupodiastoli[a-zA-Z0-9]+)_/$1\t/;
	#$line =~ s/(b\d\d\d\d)_/$1\t/; 
	#$line =~ s/(c\d\d\d\d)_/$1\t/;
	#$line =~ s/(CP\d\d\d\d)_/$1\t/;
	$line =~ s/upodiastoli/_/g;
	#$line =~ s/\t/_/;
	
	#add a tab after every gene name
	$line =~ s/(id[0-9]+_[0-9]+)_/$1\t/g;

	print $line;

}

close $duplicates_file;
