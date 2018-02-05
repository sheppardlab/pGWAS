#################################################
# GWAS summarizing of associated words script   #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#! /usr/local/bin/perl

# modified from summarize_bed.pl
use lib '/home/leonardos.mageiros/cpan/lib/perl5/site_perl/5.8.8';
use Data::Dumper;
use Getopt::Std;
use Tie::IxHash;
use strict;
use warnings;

use vars qw ($opt_f $opt_k $opt_o $opt_p $opt_a $opt_b $opt_g $opt_r $opt_v);
&getopts('f:k:p:o:a:b:g:rv');

my $usage = <<_EOH_;

## Required:
# -f assomap_v4_cf_core_clinical_results_pval/results_ST-21_s5_s12_w30.csv
# -k ST-21_s5_s12_w30 (prefix of .bed file for output)
# -p 0.0001
# -o assomap_v4_cf_core_clinical_p1E-4

# -a ~/backupsamma/SheppardLab/Campylo1/seq29.embl.cds.merged.woPAN.bed

# [-g phenotype_cat_across_original_phenotype_files/phenotype_v6.cat_ST-21_45_other9CCs.csv]

# [-r (output .overlapGenes file for reversely-associated words as well) ]

# [-v]

_EOH_
;

#
# env
#
my $bin_formatdb = "formatdb";
my $bin_blastall = "blastall";

my $bin_closestBed = "~/bin/closestBed";
my $PL_ACROSS = "summarize_assoc_words_across_CCs_freq_change.pl";

if (`ls $bin_closestBed` eq "") {
  die "Error: $bin_closestBed doesn't exist";
}

if (`which $bin_blastall` eq "") {
  die "Error: $bin_blastall doesn't exist";
}

if (`which $bin_formatdb` eq "") {
  die "Error: $bin_formatdb doesn't exist";
} 

#
# IN
#

my $inCsvFilePath  = $opt_f or die $usage;
my $key            = $opt_k or die $usage;
my $p_cutoff       = $opt_p or die $usage;

my $outDir         = $opt_o or die $usage;
   $outDir         =~ s/\/$//g;
if (! -d $outDir) {
  mkdir $outDir;
}

my $inCsvFname = $inCsvFilePath;
   $inCsvFname =~ s/^.*\///g;

my $group_name = $key;
   $group_name =~ s/results_//g;
   $group_name =~ s/_.*$//g;

my $annoBedFile = $opt_a or die $usage;
if (! -f $annoBedFile) {
  die "Error: $annoBedFile doesn't exist";
}

#my $blastCDSFile= $opt_b or die $usage;
#if (! -f $blastCDSFile) {
#  die "Error: $blastCDSFile doesn't exist";
#} else {
#  my @arr_suffix = (
#    ".nhr"
#   ,".nin"
#   ,".nsd"
#   ,".nsi"
#   ,".nsq"
#  );
#  foreach my $each_suffix (@arr_suffix) {
#    my $each_formatdb_file = $blastCDSFile . $each_suffix;
#    if (! -f $each_formatdb_file) {
#      my $cmd = "$bin_formatdb -i $blastCDSFile -n $blastCDSFile -o T -p F";
#      print($cmd);
#      system($cmd);
#      last;
#    }
#  }
#}

my $across_group_pheno_csv = "";
if ($opt_g) {
  $across_group_pheno_csv = $opt_g;

  if (! -f $across_group_pheno_csv) {
    die "Error: $across_group_pheno_csv doesn't exist";
  }
}


#
# global
#

#tie my %hash, 'Tie::IxHash';

my $cmd = "";

my %hash_word_info = (); # info in the result.csv file (e.g., presence/absence in each strain)

my $Csv_header = `grep "^word," $inCsvFilePath`;
chomp($Csv_header);

my @arr_Csv_header = split(/,/, $Csv_header);
my $header_presence_absence_each_strain = $arr_Csv_header[11];
for (my $i=12; $i<scalar(@arr_Csv_header); $i++) {
  $header_presence_absence_each_strain .= "\t" . $arr_Csv_header[$i];
}

my $header_closestBed  = "start_word";
   $header_closestBed .= "\t" . "end_word";
   $header_closestBed .= "\t" . "word_info";
   $header_closestBed .= "\t" . "positive_association";
   $header_closestBed .= "\t" . "start_gene\tend_gene";
   $header_closestBed .= "\t" . "locus";
   $header_closestBed .= "\t" . "strand";
   #$header_closestBed .= "\t" . "word_strand";
   $header_closestBed .= "\t" . "description";
   $header_closestBed .= "\t" . "distance_to_locus";
   $header_closestBed .= "\t" . "pval\tscore\tOR";
   $header_closestBed .= "\t" . $header_presence_absence_each_strain;

my $host1 = $Csv_header;      # Y=1
   $host1 =~ s/^.*host1\(//g;
   $host1 =~ s/\).*$//g;
   chomp($host1);
my $host2 = $Csv_header;      # Y=0
   $host2 =~ s/^.*host2\(//g;
   $host2 =~ s/\).*$//g;
   chomp($host2);

my %hash_pos_word_check = ();

my $out_bed = "$outDir/" . $key . ".filtered.ucsc";
open(OUT_BED_UCSC, "> $out_bed");
print OUT_BED_UCSC "track name='$key' visibility=1 itemRgb='On'\n";

my $out_bed_host1 = "$outDir/" . $key . ".filtered.$host1.bed";
open(OUT_BED_HOST1, "> $out_bed_host1");
print OUT_BED_HOST1 "track name='$key.$host1' visibility=1 itemRgb='On'\n";

my $out_bed_host2 = "$outDir/" . $key . ".filtered.$host2.bed";
open(OUT_BED_HOST2, "> $out_bed_host2");
print OUT_BED_HOST2 "track name='$key.$host2' visibility=1 itemRgb='On'\n";

my $out_csv = "$outDir/" . $inCsvFname . ".filtered";
open(OUT_CSV, "> $out_csv");
print OUT_CSV $Csv_header;

#my $out_unmapped_bed = $key . ".unmapped.filtered.ucsc";
#open(OUT_UNMAPPED_BED, "> $out_unmapped_bed");
#print OUT_UNMAPPED_BED "track name='$key' visibility=1 itemRgb='On'\n";

my $out_unmapped_bed_host1 = "$outDir/" . $key . ".unmapped.filtered.$host1.bed";
open(OUT_UNMAPPED_BED_HOST1, "> $out_unmapped_bed_host1");
print OUT_UNMAPPED_BED_HOST1 "track name='$key.$host1' visibility=1 itemRgb='On'\n";

my $out_unmapped_bed_host2 = "$outDir/" . $key . ".unmapped.filtered.$host2.bed";
open(OUT_UNMAPPED_BED_HOST2, "> $out_unmapped_bed_host2");
print OUT_UNMAPPED_BED_HOST2 "track name='$key.$host2' visibility=1 itemRgb='On'\n";

my $out_csv_unmapped = "$outDir/" . $inCsvFname . ".unmapped.filtered";
open(OUT_CSV_UNMAPPED, "> $out_csv_unmapped");
print OUT_CSV_UNMAPPED $Csv_header;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# start
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
open(IN, $inCsvFilePath);

while (my $line = <IN>) {
  if ($line =~ /^word/) {
    next;
  }
  chomp $line;
  
  #if ($line =~ /$key/) {
    #(my $stdout, my $loc, my $others_csv) = split(/\t/, $line);

    my @arr_csv = split(/,/, $line);
    # 0	word
    # 1	pvalue
    # 2	score
    # 3	OR
    # 4	whichhost
    # 5	host1(poultry)_w1[total=34]
    # 6	host2(chicken)_w1[total=35]
    # 7	host1(poultry)_w0
    # 8	host2(chicken)_w0
    # 9	location
    # 10	match
    # 11 ... presence/absence in each strain

    #
    # first, filtering by pval
    #
    my $pval  = $arr_csv[1];
    if ($p_cutoff <= $pval) {
      next;
    }
    
    my $loc   = $arr_csv[9];
    my $match = $arr_csv[10];

    #
    # prepare 
    #
    my $word = $arr_csv[0];
    my $whichhost = $arr_csv[4];

    my $word_reverse_complement = reverse $word;
       $word_reverse_complement =~ tr/ACGTacgt/TGCAtgca/;
    # skip reverse complement word because it is the same
    if (defined($hash_pos_word_check{$loc}{$word_reverse_complement}) ||
        defined($hash_pos_word_check{$loc}{$word}) ) {
      next;
    }
    $hash_pos_word_check{$loc}{$word} = 1;
    
    my $start = $loc - 14;
    my $end   = $loc + 15;
    
    my $a11 = $arr_csv[5] > 0 ? $arr_csv[5] : 1;
    my $a00 = $arr_csv[8]  > 0 ? $arr_csv[8]  : 1;
    my $a10 = $arr_csv[6] > 0 ? $arr_csv[6] : 1;
    my $a01 = $arr_csv[7]  > 0 ? $arr_csv[7]  : 1;
    my $OR    = sprintf("%0.1f",$a11 * $a00 / ($a10 * $a01));

    my $score = ($arr_csv[5]+$arr_csv[8]) - ($arr_csv[6]+$arr_csv[7]);
    
    my $id    = sprintf("%s_%s_%s_%s_%s_m%s_%s_%s_%s_comp(%s_%s)_assoc(%s)",
                        $word,$arr_csv[5],$arr_csv[6],$arr_csv[7],$arr_csv[8],$arr_csv[10],$score,$OR,$pval,$host1,$host2,$whichhost);
    my $strand = "+";
    
    my $rgb   = "";
    
    # 0,255,0   = green
    # 0,0,255   = blue
    # 255,0,255 = pink
    # 255,140,0 = orange
    # 0,0,0     = black
    # 255,0,0   = red
    if ($whichhost eq $host1) {
      $rgb = "255,0,0"; # 
    } elsif ($whichhost eq $host2) {
      $rgb = "0,255,0"; # green
    } else {
      # white
      $rgb = "255,255,255";
    }

    my $abs_score = abs($score);

    #
    # mapped
    #
    if ($loc =~ /^[-0-9]+/ && $loc > -1) {
      #
      # for perfect match words
      #
      if ($match == 0) {
        my $line_presence_absence_each_strain = $arr_csv[11];
        for (my $i=12; $i<scalar(@arr_csv); $i++) {
          $line_presence_absence_each_strain .= "\t" . $arr_csv[$i];
        }
        $hash_word_info{$word}{'presence_absence_each_strain'} = $line_presence_absence_each_strain;
        $hash_word_info{$word}{'pval'} = $pval;
        $hash_word_info{$word}{'score'} = $score;
        $hash_word_info{$word}{'OR'} = $OR;
        
        #
        # output
        #
        print OUT_BED_UCSC "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
        
        if ($whichhost eq $host1) {
          print OUT_BED_HOST1 "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
        }
        if ($whichhost eq $host2) {
          print OUT_BED_HOST2 "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
        }
        
        print OUT_CSV $line . "\n";
      }
    #
    # unmapped
    #
    } else {
      # chr 65773 65866 50  0 . 65773 65866 R,G,B
      #print OUT_UNMAPPED_BED         "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
      if ($whichhost eq $host1) {
        print OUT_UNMAPPED_BED_HOST1 "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
      }
      if ($whichhost eq $host2) {
        print OUT_UNMAPPED_BED_HOST2 "chr\t$start\t$end\t$id\t$abs_score\t$strand\t$start\t$end\t$rgb\n";
      }
      print OUT_CSV_UNMAPPED $line . "\n";
    }

  #}

}
close(IN);
close(OUT_BED_UCSC);

close(OUT_BED_HOST1);
close(OUT_BED_HOST2);

close(OUT_UNMAPPED_BED_HOST1);
close(OUT_UNMAPPED_BED_HOST2);

close(OUT_CSV);
close(OUT_CSV_UNMAPPED);


# - - - - - - - - - - - - - - - - - - - - 
# closestBed as pipe
# - - - - - - - - - - - - - - - - - - - - 

# 0	chr
# 1	5095
# 2	5124
# 3	ATGCATATTGGTAATGAGATTGGATCCATG_15_9_2_9_m0_13_7.5_6.2e-05_comp(s5_s4)
# 4	(s5)
# 5	13
# 6	+
# 7	5095
# 8	5124
# 9	255,0,0
# 10	chr
# 11	4916
# 12	5257
# 13	CAMP0004_Cj0004c_Cj0004c
# 14	0
# 15	-
# 16	4916
# 17	5257
# 18	0,0,0
# 19	putative periplasmic protein
# 20	0

my @arr_target_host_assoc_bed = ();
push(@arr_target_host_assoc_bed, $out_bed_host1);
if ($opt_r) {
  push(@arr_target_host_assoc_bed, $out_bed_host2);
}

foreach my $each_target_bed (@arr_target_host_assoc_bed) {

  # final output of interest
  my $out_closestBed_host_assocGenes     = $each_target_bed . ".overlapGenes";
  my $out_closestBed_host_assocGenes_fas = $each_target_bed . ".overlapGenes.fas";

  my %hash_host_assocGenes = ();

  open(OUT_CLOSESTBED_HOST_ASSOCGENES, "> $out_closestBed_host_assocGenes");
  print OUT_CLOSESTBED_HOST_ASSOCGENES $header_closestBed . "\n";

  $cmd = "$bin_closestBed -a $each_target_bed -b $annoBedFile -d | sort -k 2,2 -u -n |";
  open(CLOSESTBED, "$cmd");
  while (my $line = <CLOSESTBED>) {
    chomp $line;
    $line =~ s/_assoc/\t/g;
    my @arr_line = split(/\t/, $line);

    my $word_info = "";
    my $word_seq = $arr_line[3];
       $word_seq =~ s/_.*$//g;
    my $word_assoc_info = $arr_line[3];
       $word_assoc_info =~ s/(^[ATGCatgz]+_)//g;

#    #
#    # 24/10/2013
#    # get strand of words, 
#    # and output word sequences in the same direction as genes on which they locate
#    # 
#    my $word_strand = "";
#    
#    my $cmd_blast = "echo '$word_seq' | $bin_blastall -p blastn -d '$blastCDSFile' -m 8 2> /dev/null | head -1";
#    if ($opt_v) {
#      print $cmd_blast . "\n";
#    }
#    my $tophit_blast = `$cmd_blast`;
#    chomp($tophit_blast);
#    
#    if ($tophit_blast ne "") {
#      my @arr_tophit_blast = split(/\t/, $tophit_blast);
#      #
#      # tabular format
#      #   http://www.pangloss.com/wiki/Blast
#      # 
#      my $hit_start = $arr_tophit_blast[8];
#      my $hit_end   = $arr_tophit_blast[9];
#      if ($hit_start > $hit_end) {
#        $word_seq = reverse $word_seq;
#        $word_seq =~ tr/ACGTacgt/TGCAtgca/;
#      }
#      $word_strand = "same_strand_as_gene";
#    } else {
#      $word_strand = "unknown_strand (no hit on CDS by blast)";
#    }

    $word_info = $word_seq . "_" . $word_assoc_info;

    my $out_line  =        $arr_line[1];   # start_word
       $out_line .= "\t" . $arr_line[2];   # end_word
       $out_line .= "\t" . $word_info;     # word_info
       $out_line .= "\t" . $arr_line[4];   # positive_association
       $out_line .= "\t" . $arr_line[11];  # start_gene
       $out_line .= "\t" . $arr_line[12];  # end_gene
       $out_line .= "\t" . $arr_line[13];  # locus
       $out_line .= "\t" . $arr_line[15];  # strand
       #$out_line .= "\t" . $word_strand; 
       $out_line .= "\t" . $arr_line[19];  # description
       $out_line .= "\t" . $arr_line[20];  # distance_to_locus

    $hash_host_assocGenes{$arr_line[13]} = 1;

    my $word = $arr_line[3];
       $word =~ s/_.*$//g;

    if (!defined($hash_word_info{$word}{'pval'})) {
      print Dumper(\@arr_line);
      die "Error: $word is not defined in hash_word_info";
    }

    print OUT_CLOSESTBED_HOST_ASSOCGENES $out_line;
    print OUT_CLOSESTBED_HOST_ASSOCGENES "\t" . $hash_word_info{$word}{'pval'};
    print OUT_CLOSESTBED_HOST_ASSOCGENES "\t" . $hash_word_info{$word}{'score'};
    print OUT_CLOSESTBED_HOST_ASSOCGENES "\t" . $hash_word_info{$word}{'OR'};
    print OUT_CLOSESTBED_HOST_ASSOCGENES "\t" . $hash_word_info{$word}{'presence_absence_each_strain'};
    print OUT_CLOSESTBED_HOST_ASSOCGENES "\n";

  }
  close(CLOSESTBED);

  my @arr_host1_assocGenes = keys %hash_host_assocGenes;
  print "$out_closestBed_host_assocGenes\t" . scalar(@arr_host1_assocGenes) . " genes\n";
  
  close(OUT_CLOSESTBED_HOST_ASSOCGENES);

}


# - - - - - - - - - - - - - - - - - - - - 
# across CCs
#   output to $out_bed_host1.across
# - - - - - - - - - - - - - - - - - - - - 

if ($across_group_pheno_csv ne "") { # $opt_g

  if (! -f $PL_ACROSS) {
    die "Error: $PL_ACROSS doesn't exist";
  }

  $cmd = "qsub -cwd -l sjob <<< 'perl $PL_ACROSS -b $out_bed_host1 -p $across_group_pheno_csv -g $group_name' ";
  print($cmd);
  system($cmd);

}



#
# TODO
#   select representative associated genetic variation from the list of words
#   give a file for iTOL (although it can be done manually)
#
