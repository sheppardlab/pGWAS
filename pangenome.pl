#################################################
# Pangenome list creation script                #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/bin/perl

use lib '~/lib/site_perl/5.14.2';
use lib '/home/leonardos.mageiros/lib/site_perl/5.14.2';
use Data::Dumper;
use Getopt::Std;
use Tie::IxHash;
use strict;
use warnings;

use vars qw ($opt_r $opt_c $opt_t $opt_l $opt_o $opt_v);
&getopts('r:c:t:l:ov');

my $usage = <<_EOH_;

## Required:
# -r Ecoli_K-12_MG1655_chromosome.txt
# -c Ecoli_ (prefix of *.txt files)
# -t 70
# -l 50

# [-o]

# [-v]

_EOH_
;

#
# IN
#

my $refFile   = $opt_r or die $usage;
my $compPrefix = $opt_c or die $usage;
my $threshold_iden = $opt_t or die $usage;
my $threshold_len = $opt_l or die $usage;

my $opposit_flag = 0;
if ($opt_o) {
  $opposit_flag = 1;
}

system("perl -i -pe 's/\r//g' $refFile");

#if ($opt_v) {
#
#}

# 0     id
# 1     data_type
# 2     allele_id_format
# 3     description
# 4     length
# 5     length_varies
# 6     coding_sequence
# 7     flag_table
# 8     main_display
# 9     isolate_display
# 10    query_field
# 11    analysis
# 12    reference_sequence

tie my %out_hash_union, 'Tie::IxHash';
#tie my %hash_headerIndex2Name, 'Tie::IxHash';

my $out_unionTxt = "union.txt";
system("cat $refFile > $out_unionTxt\n");
open(OUT_UNION_TXT, ">> $out_unionTxt");

my $out_unionFas = "union.fas";
open(OUT_UNION_FAS, "> $out_unionFas");

my $out_union_log = "union.log";
open(OUT_UNION_LOG, "> $out_union_log");

my $out_duplicates_log = "duplicates.log";
open(OUT_DUPLICATES_LOG, "> $out_duplicates_log");

my $log_header = "added_id\t%aln_length\t%identity\ttop_hit_id\n";
print OUT_UNION_LOG       $log_header;
print OUT_DUPLICATES_LOG  $log_header;

my %hash_id2sub = ();

#
# read reference info
#
my %hash_refseq = ();
open(REF, $refFile);

my $line_header = <REF>;
#my @arr_line_header = split(/\t/, $line_header);
#for (my $i=0; $i<scalar(@arr_line_header); $i++) {
#  print $hash_headerIndex2Name{ $i } = $arr_line_header[$i];
#}

while (my $line = <REF>) {
  chomp $line;
  my @arr_line = split(/\t/, $line);

  $arr_line[0] =~ s/^([0-9a-zA-Z-_]+).*/$1/g;

  my $id = $refFile . "_" . $arr_line[0] . "_" . $arr_line[3];
  my $seq = $arr_line[12];

  $id =~ s/ /_/g; # remove space

  if (defined($hash_id2sub{$id})) {
    $hash_id2sub{$id}++;
    $id .= "_" . $hash_id2sub{$id};
  }

  print OUT_UNION_FAS ">$id\n";
  print OUT_UNION_FAS "$seq\n";

  # record this finished id
  if (!defined($hash_id2sub{$id})) {
    $hash_id2sub{$id} = 1;
  }
}
close(REF);

close(OUT_UNION_FAS);

my $opt_formatdb = "-o F -p F";
my $cmd_formatdb = "formatdb -i $out_unionFas -n $out_unionFas $opt_formatdb\n";
print($cmd_formatdb);
system($cmd_formatdb);

#
# start comparison loop
#
my @arr_compFiles = ();
if (-f $compPrefix) {
  push(@arr_compFiles, $compPrefix);
} else {
  @arr_compFiles = glob("$compPrefix*.txt");
}

foreach my $each_compFile (@arr_compFiles) {
  if ($each_compFile ne $refFile) {

    my $i=0;
    open(COMP, $each_compFile);
    while (my $line = <COMP>) {
      if ($line =~ /^id\t/) {
        next;
      }
      chomp $line;
      my @arr_line = split(/\t/, $line);

      if ($i % 100 == 0) {
        print "Examining $i th entry ... \n";
      }

      $arr_line[0] =~ s/^([0-9a-zA-Z-_]+).*/$1/g;

      my $id = $each_compFile . "_" . $arr_line[0] . "_" . $arr_line[3];
      my $query_seq = $arr_line[12];
      my $query_len = length($query_seq);

      my $cmd_blast = "echo '" . $query_seq . "' | blastall -p blastn -d '$out_unionFas' -m 8  | head -1 2>/dev/null";
      if ($opt_v) {
        #print "$cmd_blast\n";
      }
      my $tophit_blast_ref = `$cmd_blast`;

      if ($tophit_blast_ref ne "") {
        my @arr_tophit_blast_ref = split(/\t/, $tophit_blast_ref);
        #
        # tabular format
        #   http://www.pangloss.com/wiki/Blast
        #
        #   0	Query
        #   1	Subject
        #   2	% id
        #   3	alignment length (including gaps)
        #   4	mistmatches
        #   5	gap openings
        #   6	q.start
        #   7	q.end
        #   8	s.start
        #   9	s.end
        #   10	e-value
        #   11	bit score

        #
        my $top_hit_id = $arr_tophit_blast_ref[1];
        my $per_iden  = $arr_tophit_blast_ref[2];
        
        my $aligned_len  = $arr_tophit_blast_ref[3]; 
        my $per_aligned_len = sprintf("%0.2f",($aligned_len/$query_len)*100);
        
        if ($opt_v) {
          print "$id\t$per_iden\t$per_aligned_len\t$top_hit_id\n";
        }

        # new locus which doesn't exist in the reference
        my $new_flag = 0;
        if ($per_aligned_len < $threshold_len) {
          $new_flag = 1;
        } elsif ($per_iden < $threshold_iden) {
          $new_flag = 1;
        }
        
        if ( ($new_flag == 1 && $opposit_flag == 0) || ($new_flag == 0 && $opposit_flag == 1) ) {
          print "New locus: $id (%aln_length=$per_aligned_len ($aligned_len/$query_len), %identity=$per_iden, top_hit_id=$top_hit_id)\n";
          print OUT_UNION_LOG "$id\t$per_aligned_len\t$per_iden\t$top_hit_id\n";

          $id =~ s/ /_/g; # remove space

          if (defined($hash_id2sub{$id})) {
            $hash_id2sub{$id}++;
            $id .= "_" . $hash_id2sub{$id};
          }

          # add it to the union fas
          open(OUT_UNION_FAS, ">> $out_unionFas");
          print OUT_UNION_FAS ">$id\n";
          print OUT_UNION_FAS "$query_seq\n";
          close(OUT_UNION_FAS);
          
          # add it to the union txt
          print OUT_UNION_TXT $line . "\n";

          # update the union blastDB
          my $cmd_formatdb = "formatdb -i $out_unionFas -n $out_unionFas $opt_formatdb\n";
          print($cmd_formatdb);
          system($cmd_formatdb);
          
          # record this finished id
          if (!defined($hash_id2sub{$id})) {
            $hash_id2sub{$id} = 1;
          }
        } else {
          print OUT_DUPLICATES_LOG "$id\t$per_aligned_len\t$per_iden\t$top_hit_id\n";
        }
        
      } # if ($tophit_blast_ref ne "") {

      $i++;
    }
    close(COMP);

  }
}

close(OUT_UNION_TXT);
close(OUT_UNION_LOG);
close(OUT_DUPLICATES_LOG);