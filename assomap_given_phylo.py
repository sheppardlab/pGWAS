#################################################
# GWAS execution script                         #
# Script from the Sheppard Lab pGWAS pipeline   #
# Mageiros, Meric et al. July 2016              #
# For more information: s.k.sheppard@bath.ac.uk #
# http://www.sheppardlab.com/                   #
#################################################

#!/usr/local/bin/python -u
# coding: EUC-JP

import re
import numpy
import getopt
import time
import sys;
import pprint;
import glob;
import csv;
import Bio.SeqIO;
import Bio.Phylo;
import shelve;
from math import exp
import random
import collections

import rlcompleter, readline
readline.parse_and_bind("tab: complete")

pp = pprint.PrettyPrinter(indent=4)
#signif_cutoff = 0.0001

def usage():
	print "usage: ./assomap_given_phylo.py [OPTION]"
	print "[OPTION]"

	print " -s 0.0001"
	print " -d Campylo2_contigs/ST-45 (where *.ctg are stored)"
	print " -f 30"
	print " -p phenotype_tree_files/phenotype_5stages_ST-45.csv (int_id,prefix_of_ctg_file,host_str,any...)" 
	print " -t phenotype_tree_files/ST-45.int.nwk (OUT name must be replaced with integer id (0,...,n-1) in the phenotype file)"
	print " -r ../Campylo1/contigs/seq29.fas"
	print " -c s45_s123 (1,0)"
	print " -d Ecoli_ECOR_GMB_contigs (where *.ctg are stored)"

#	print " -f 30"
#	print " -p Ecoli_ECOR_GMB_contigs/171S_SNPs.Base.int.phenotype.csv"
#	print " -t Ecoli_ECOR_GMB_contigs/171S_SNPs.Base.int.nwk (OUT name must be replaced with integer id (0,...,n-1) in the phenotype file)"
#	print " -r Ecoli_ECOR_GMB_contigs/K12_NC_000913.fas"
#	print " -c 1,2"
	sys.exit()

def approxkeys(key,level):
  if level==0:
    return [key];
  feat=len(key);
  keys1=approxkeys(key,level-1);
  keys=[];
  for key in keys1:
    for j in range(0,feat):
      for letter in ['A','C','G','T']:
        keys.append(key[0:j]+letter+key[j+1:]); # note that 0:0 means nothing
  return keys;

#############################################################################
score_range = 400

#############################################################################
try:
	shortopt = "s:d:f:p:t:r:c:"
	longopt = ["signif_cutoff=","dirname=", "featlength=", "phenofile=", "treefile=", "ref_fastafile=", "compare="]
	opts, args = getopt.getopt(sys.argv[1:], shortopt, longopt)
except getopt.GetoptError:
	usage()

signif_cutoff = None
dirname = None
feat = None
phenofile = None
treefile = None
ref_fastafile = None
compare = None

verbose = False

for o, a in opts:
	if o in ("-s", "--signif_cutoff"):
		signif_cutoff = a
	elif o in ("-d", "--dirname"):
		dirname = a
	elif o in ("-f", "--featlength"):
		feat = a #length of a word
	elif o in ('-p', '--phenofile'):
		phenofile = a
	elif o in ('-t', '--treefile'):
		treefile = a
	elif o in ('-r', '--ref_fastafile'):
		ref_fastafile = a
	elif o in ('-c', '--compare'):
		compare = a

if signif_cutoff == None or dirname == None or feat == None or phenofile == None or treefile == None or ref_fastafile == None or compare == None :
	usage()

signif_cutoff = float(signif_cutoff)
print(signif_cutoff)

ST = re.sub("^.*/",'',dirname)

feat = int(feat)

compare_binary_str_list = compare.split('_')
# compare_binary_str_list[0] <=> host '1' (previously)
# compare_binary_str_list[1] <=> host '2' (previously, better to write '0')

#binary_str2int_dict = {};
#binary_str2int_dict[compare_binary_str_list[0]] = 1
#binary_str2int_dict[compare_binary_str_list[1]] = 0

#############################################################################
print('Reading the host file')
d_intId=[];
d=[];
#
# d_intId[i]              : intId of each sample used in the tree 
#                           (previously, int(d[i][len(dirname)+1:-4]))
#
# d[i]                    : path of each .ctg file 
# d[i][len(dirname)+1:-4] : prefix of .ctg file (BIGS ID)
#

#
# phenotype.csv
#   intId used in the tree         => d_intId[]
#   prefix of .ctg file (BIGS ID)  => d[]
#   host_str
#   any string
#
# you can use BIGS ID for both d_intId[] and d[] (confirmed in July, 2013)
# by specifying BIGS ID for both 1st and 2nd columns, and for OUT in the tree
#
# you also can use ClonalFrame ID for d_intId[],
#
handle=open(phenofile) 
f=csv.reader(handle) 

sampleIntId2host_dict={};
for record in f:
  d.append(dirname+'/'+record[1]+'.ctg')
  d_intId.append(int(record[0]))
  sampleIntId2host_dict[int(record[0])]=record[2];

handle.close()
n=len(d); #Number of isolates

#############################################################################
print('Monte-Carlo test')
sim_data_dict={}#simulated dataset 
repMC=1000000#Number of replicates

tree=Bio.Phylo.read(treefile,'newick') 

allclades = list(tree.find_clades(order='level'))
sumbralen=0.0 # sum of branch length
for clade in allclades:
  sumbralen+=clade.branch_length

rate=1.0/sumbralen # mut rate (equal on each branch(=clade))

# add attribute seq (random 0/1 pattern) to root () in the tree
tree.root.seq=[0]*repMC # = allclades[0].seq
for i in range(0,repMC):
  tree.root.seq[i]=int(round(random.random())) # round => 0 or 1

#
# prepare 
#   sim_data_dict[] 
#        (sample_id => array of 0/1 whose size is repMC)
#         s.t. samples close to each other in phylogeny will share more 0/1 pattern
#
for clade in allclades: # clade = branch
  # start from root = allclades[0]
  # >>> allclades[0]
  # Clade(branch_length=1.0)
  #
  for child in clade:
    # >>> len(allclades[0])
    # 3
    # >>> allclades[0][0]
    # Clade(branch_length=0.00012377, confidence=0.78) branch to the lower part (see the tree in MEGA)
    # >>> allclades[0][1]
    # Clade(branch_length=0.00217106, confidence=33.0) OUT
    # >>> allclades[0][2]
    # Clade(branch_length=0.00636023, confidence=11.0) OUT
    # 
    # >>> len(allclades[1])
    # 2
    # >>> allclades[1][0]
    # Clade(branch_length=0.00053612, confidence=0.97)
    # >>> allclades[1][1]
    # Clade(branch_length=0.00114403, confidence=1.0)
    #
    #
    #
    dist=child.branch_length
    psame=(1.0+exp(-2.0*dist*rate))/2.0 # prob to keep the same pattern from parent to child
    # 
    #   2.0*dist*rate : param of Poisson (number of mut (both ways))
    #   
    #   Question: words were simulated to evolve through a process of gain and loss along the branches of the phylogeny.
    #             reference
    #
    
    # prepare child.seq (child = evolved) array evolved from the parent clade 
    child.seq=[0]*repMC 
    for i in range(0,repMC):
      if random.random()<psame:child.seq[i]=clade.seq[i] # keep the presence/absence pattern of parent clade
      else:child.seq[i]=1-clade.seq[i]                   # flip the presence/absence pattern (gain or loss)
    
    # set child.seq array to sim_data_dict[sample_id]
    if child.is_terminal():sim_data_dict[int(child.confidence)]=child.seq
                                           # child.confidence  =sample_id (OUT name of integer id) of the child branch
#
#
#
print('Convert host and sim_data from dict to array ordered as the sequence files')
host=[];
sim_data=[];
for i in range(0,n):
  id=d_intId[i] # i = id only if intId is defined as 0,...,n-1
  host.append(sampleIntId2host_dict[id])
  sim_data.append(sim_data_dict[id])

#
#
#
print('Compute typical scores from simulation')
subset=[] 
for i in range(0,n):
  if (host[i]==compare_binary_str_list[0]) | (host[i]==compare_binary_str_list[1]):subset.append(i); 
  # record samples of host 1 and host 2
  # i.e. infomraiton of host (phenotype) of each sample is fixed on the phylongey
  #      presence/absence pattern of a word among the samples is varied by simulation

scores=[0]*repMC
scoresOR=[0]*repMC
for r in range(0,repMC): # for the samples of host 1 and host 2,
  score=0
  n1=0;n2=0;n3=0;n4=0; # only for OR and its null distribution
  #          word==1  word!=1
  # host==1    n1      n2
  # host!=1    n3      n4
  #
  for i in subset: # count assocation of presence/absence (conceptually, 2x2)
    if ((host[i]==compare_binary_str_list[0]) & (sim_data[i][r]==1)) | ((host[i]!=compare_binary_str_list[0]) & (sim_data[i][r]==0)):score=score+1;
    else:score=score-1

    if   ((host[i]==compare_binary_str_list[0]) & (sim_data[i][r]==1)):
      n1=n1+1;
    elif ((host[i]==compare_binary_str_list[0]) & (sim_data[i][r]==0)):
      n2=n2+1;
    elif ((host[i]!=compare_binary_str_list[0]) & (sim_data[i][r]==1)):
      n3=n3+1;
    elif ((host[i]!=compare_binary_str_list[0]) & (sim_data[i][r]==0)):
      n4=n4+1;

  if (n1==0):n1=1;
  if (n2==0):n2=1;
  if (n3==0):n3=1;
  if (n4==0):n4=1;

  scores[r]=score
  scoresOR[r]=float( (n1*n4) )/(n2*n3) 

#
# null distribution of score was obtained above
#

#f=open('scores.csv','w')
#for i in range(0,repMC):f.write(str(scores[i])+'\n');
#f.close()
print u"scores: average=%0.1f, std=%0.1f, ptp=%s, min=%s, max=%s" % ( numpy.average(scores), numpy.std(scores), str(numpy.ptp(scores)), str(numpy.min(scores)), str(numpy.max(scores)) )
# average: around 0
# range: small pop structure => small std of null dist of scores


#############################################################################
print('Create hash table')
ht = {};
ht_each_host = collections.defaultdict(lambda:collections.defaultdict(int));
for i in range(0,n):
  print 'Adding file '+d[i]+' to hash table...'
  handle=open(d[i], 'r')
  f=Bio.SeqIO.parse(handle, 'fasta')  
  thisone='0'*i+'1'+'0'*(n-i-1)
  thisone_csv='0,'*i+'1'+',0'*(n-i-1)
  j=0;
  for record in f:
    j=j+1;
    #if j>1:break; # must be removed
    l=len(record.seq)
    seq = str(record.seq)+str(record.seq.reverse_complement())
    seq = seq.upper() # added by Koji

    #for k in range(0,l-feat)+range(l,2*l-feat): # range(a,b) does not include the b.  Currently, a word at the end of a contig is not included in this hash table.
    for k in range(0,l-feat+1)+range(l,2*l-feat+1): # modified by Koji
      word=seq[k:k+feat]

      #ht_each_host[word][d[i][len(dirname)+1:-4]] = 1 # only for debugging

      # ht is shared among the all samples to record words
      #   word => n-digit 0/1 which expresses who has the word
      if not(ht.has_key(word)):
        ht[word]=thisone;
      elif ht[word][i]=='0':
        ht[word]=ht[word][0:i]+'1'+ht[word][i+1:]
f.close()

#
#
#
print 'Converting hash table to a 2D array...'
all=ht.items();
# [('CAAAAGCTTCAAAAAGTTCAAGGGTATTCT', '000000001000000000000010000000000000000000000000001'), ('CTTGATCGCTTAGATGAAGAACCGACTCGA', '001000000000000000000000000000000000000000000000000'), ('TAGAAAATTTAGATTTAAAGACTTTAAATT', '000100000000000000110000000000000000001000000000000'), ('ATTTCTACGCCTTTATCTTTAGCACGCTCT', '000000000000000000000000000000000000100000100100000'), ('CGCAAATTTTGTTCAAGCATAGTTTTGTAT', '010000000000000000000000100000010001010101000000000')] ...


#############################################################################
print 'Looking for location in reference'
loc={};
handle=open(ref_fastafile,'r'); 
f=Bio.SeqIO.parse(handle,'fasta')
record=f.next()#Note that only the first contig of the reference is used
l=len(record.seq)
seq=str(record.seq)+str(record.seq.reverse_complement())

for k in range(0,l-feat)+range(l,2*l-feat): # k: start of a word
  word=seq[k:k+feat]
  if not(loc.has_key(word)):
    lo=k+feat/2;       # lo: middle pos of a word, from 5 to 3 prime
    if k>=l:
      lo=2*l-(k+feat/2); # lo: middle pos of a word, from 3 to 5 prime
    loc[word]=lo;

f.close()

#############################################################################
print 'Preparing p-value of each score in the null distribution ...'

#
# then calculate p-value of observation of each value of score (from bottom score_range/2 to top score_range/2 points)
#
pval =[0.0]*score_range; # p-value of score_range/2 points for each tail
pval2=[0.0]*score_range; # p-value of score_range/2 points for each tail
for i in range(0,repMC):
  for j in range(scores[i]+score_range/2,score_range):pval [j]+=1.0/repMC # count the number of simulated scores to calculate p-value on the left (scores near 0 will be 0)
  for j in range(0,scores[i]+score_range/2+1):pval2[j]+=1.0/repMC # count the number of simulated scores to calculate p-value on the right (scores near score_range will be 0)
for s in range(0,score_range):pval[s]=min(pval[s],pval2[s])/2.0;

# e.g.,
# >>> pval[0:99] # left (these values are variable, so check pval.csv)
# [0.0, 0.0, 0.0, 0.0, 0.0, 5e-07 (p-value of score=-95, observed only once) , 5e-07, 3.749999999999994e-05, 3.749999999999994e-05, 7.049999999999986e-05, 7.049999999999986e-05, 0.00020700000000000178, 0.00020700000000000178, 0.0005905, 0.0005905, 0.0006744999999999929, 0.0006744999999999929, 0.0009499999999999697, 0.0009499999999999697, 0.0013180000000000127, 0.0013180000000000127, 0.0014440000000000295, 0.0014440000000000295, 0.0015620000000000451, 0.0015620000000000451, 0.0015845000000000481, 0.0015845000000000481, 0.001614000000000052, 0.001614000000000052, 0.0017240000000000666, 0.0017240000000000666, 0.0017555000000000708, 0.0017555000000000708, 0.0018485000000000831, 0.0018485000000000831, 0.0018880000000000884, 0.0018880000000000884, 0.0019770000000001, 0.0019770000000001, 0.0020000000000001033, 0.0020000000000001033, 0.0020095000000001045, 0.0020095000000001045, 0.0020400000000001086, 0.0020400000000001086, 0.0021655000000001252, 0.0021655000000001252, 0.002208000000000131, 0.002208000000000131, 0.0022265000000001333, 0.0022265000000001333, 0.002253500000000137, 0.002253500000000137, 0.0022860000000001412, 0.0022860000000001412, 0.0024075000000001573, 0.0024075000000001573, 0.0025515000000001765, 0.0025515000000001765, 0.0025825000000001806, 0.0025825000000001806, 0.002607500000000184, 0.002607500000000184, 0.00265350000000019, 0.00265350000000019, 0.002706500000000197, 0.002706500000000197, 0.0027255000000001996, 0.0027255000000001996, 0.0028230000000002125, 0.0028230000000002125, 0.0028745000000002193, 0.0028745000000002193, 0.0030480000000002424, 0.0030480000000002424, 0.0033860000000002872, 0.0033860000000002872, 0.0056684999999990614, 0.0056684999999990614, 0.03215050000002182, 0.03215050000002182, 0.03536650000002504, 0.03536650000002504, 0.03881100000002848, 0.03881100000002848, 0.042693000000032365, 0.042693000000032365, 0.05733450000004701, 0.05733450000004701, 0.21820249999762098, 0.21820249999762098, 0.2349759999971722, 0.2349759999971722, 0.2443084999969225, 0.2443084999969225, 0.24801899999682323, 0.24801899999682323, 0.2495264999967829, 0.2495264999967829]

# >>> pval[score_range/2:199] # right (these values are variable, so check pval.csv)
# [0.2499824999967707, 0.2500174999967707, 0.24953999999678253, 0.24953999999678253, 0.2479934999968239, 0.2479934999968239, 0.24433399999692182, 0.24433399999692182, 0.23522849999716544, 0.23522849999716544, 0.2184854999976134, 0.2184854999976134, 0.05781950000004749, 0.05781950000004749, 0.04306100000003273, 0.04306100000003273, 0.03914150000002881, 0.03914150000002881, 0.03558600000002526, 0.03558600000002526, 0.032293000000021964, 0.032293000000021964, 0.005564499999999138, 0.005564499999999138, 0.00333200000000028, 0.00333200000000028, 0.002994000000000235, 0.002994000000000235, 0.002833500000000214, 0.002833500000000214, 0.0027760000000002063, 0.0027760000000002063, 0.0026870000000001944, 0.0026870000000001944, 0.0026720000000001925, 0.0026720000000001925, 0.0026130000000001846, 0.0026130000000001846, 0.002563500000000178, 0.002563500000000178, 0.0025430000000001753, 0.0025430000000001753, 0.002509000000000171, 0.002509000000000171, 0.002375000000000153, 0.002375000000000153, 0.0022525000000001368, 0.0022525000000001368, 0.002231000000000134, 0.002231000000000134, 0.00220250000000013, 0.00220250000000013, 0.002193500000000129, 0.002193500000000129, 0.0021545000000001238, 0.0021545000000001238, 0.0020340000000001078, 0.0020340000000001078, 0.0020000000000001033, 0.0020000000000001033, 0.001991000000000102, 0.001991000000000102, 0.0019725000000000996, 0.0019725000000000996, 0.0018835000000000878, 0.0018835000000000878, 0.0018430000000000824, 0.0018430000000000824, 0.001735000000000068, 0.001735000000000068, 0.001704000000000064, 0.001704000000000064, 0.001605500000000051, 0.001605500000000051, 0.0015715000000000464, 0.0015715000000000464, 0.0015440000000000427, 0.0015440000000000427, 0.0014285000000000274, 0.0014285000000000274, 0.0012915000000000092, 0.0012915000000000092, 0.0009284999999999716, 0.0009284999999999716, 0.0006864999999999919, 0.0006864999999999919, 0.0006204999999999975, 0.0006204999999999975, 0.00023500000000000246, 0.00023500000000000246, 8.499999999999982e-05, 8.499999999999982e-05, 3.9999999999999936e-05, 3.9999999999999936e-05 (p-value of score=95, observed twice) , 1e-06, 1e-06, 0.0, 0.0, 0.0

pval_file = "pval_of_scores_%s_%s_w%s.csv" % (ST, compare, feat)
f=open(pval_file,'w') 
for i in range(-score_range/2,score_range/2):f.write(str(i)+','+str(pval[i+score_range/2])+'\n'); # p-values of scores from -score_range/2 to score_range/2
f.close()

###
scoresOR.sort(reverse=False)
cutoff_small_scoreOR = scoresOR[ int(repMC*signif_cutoff/2) ]

scoresOR.sort(reverse=True)
cutoff_large_scoreOR = scoresOR[ int(repMC*signif_cutoff/2) ]

print u"scoresOR: average=%0.1f, std=%0.1f, ptp=%s, min=%s, max=%s, cutoff_small_scoreOR=%s, cutoff_large_scoreOR=%s" % ( numpy.average(scoresOR), numpy.std(scoresOR), str(numpy.ptp(scoresOR)), str(numpy.min(scoresOR)), str(numpy.max(scoresOR)), str(cutoff_small_scoreOR), str(cutoff_large_scoreOR) )

scoresOR_file = "scoresOR_%s_%s_w%s.csv" % (ST, compare, feat)
f=open(scoresOR_file,'w')
for i in range(0,len(scoresOR)-1):f.write(str(i)+','+str(scoresOR[i])+'\n'); 
f.close()



#############################################################################
print 'Screening of significant words ...'

c2=0.0;
c4=0.0;
for i in subset:
  if (host[i]==compare_binary_str_list[0]): 
    c2=c2+1; # host==1
  else:
    c4=c4+1; # host!=1

min_p = 1
result_file = "results_%s_%s_w%s.csv" % (ST, compare, feat)
f1=open(result_file,'w')

resultOR_file = "results_OR_%s_%s_w%s.csv" % (ST, compare, feat)
#fOR=open(resultOR_file,'w')

#header='word,pvalue,score,whichhost,host1[total='+str(int(c2))+'],host2[total='+str(int(c4))+'],location,match,';
header= u"word,pvalue,score,OR,whichhost,host1(%s)_w1[total=%s],host2(%s)_w1[total=%s],host1(%s)_w0,host2(%s)_w0,location,match," % (compare_binary_str_list[0], str(int(c2)), compare_binary_str_list[1], str(int(c4)), compare_binary_str_list[0], compare_binary_str_list[1]);

for i in range(0,n):header=header+d[i][len(dirname)+1:-4]+'['+host[i]+'],' # prepare header str by using d[i] (originally, names of contigs = BIGS ID are used for the header)
header=header[0:-1]
f1.write(header+'\n')
#fOR.write(header+'\n')

for a in range(0,len(all)):
  # word => n-digit 0/1 which expresses who has the word
  key=all[a][0]
  pat=all[a][1]

  patstr=''
  for i in range(0,n):patstr=patstr+pat[i]+',' # can be subset
  patstr=patstr[0:-1]

  c1=0.0;c3=0.0;
  a2=0.0;a4=0.0;

  score=0
  for i in subset:
    if ((host[i]==compare_binary_str_list[0]) & (pat[i]=='1')):c1=c1+1; # host==1 && word==1
    if ((host[i]!=compare_binary_str_list[0]) & (pat[i]=='1')):c3=c3+1; # host!=1 && word==1

    if ((host[i]==compare_binary_str_list[0]) & (pat[i]!='1')):a2=a2+1; # host==1 && word!=1
    if ((host[i]!=compare_binary_str_list[0]) & (pat[i]!='1')):a4=a4+1; # host!=1 && word!=1

    #
    #          word==1  word!=1
    # host==1    c1(0)   19       c2(19)
    # host!=1    c3(2)    3       c4(5)
    #
    # score(-18) = (0+3)-(19+2)
    #

    # count (host==1 & word==1) or (host!= 1 & word!=1) i.e. association between presence/absence of the word and host1/host2
    if ((host[i]==compare_binary_str_list[0]) & (pat[i]=='1')) | ((host[i]!=compare_binary_str_list[0]) & (pat[i]=='0')):score=score+1;
    # the opposite i.e. association between presence/absence of the word and host2/host1
    else:score=score-1;

  c1_nz=c1;c3_nz=c3;
  a2_nz=a2;a4_nz=a4;

  if (c1_nz==0):c1_nz=1;
  if (a2_nz==0):a2_nz=1;
  if (c3_nz==0):c3_nz=1;
  if (a4_nz==0):a4_nz=1;
  scoreOR = float( (c1_nz*a4_nz) )/(a2_nz*c3_nz)

  if pval[score+score_range/2]<min_p: # note: score can be negative
    min_p = pval[score+score_range/2]
    print u"lowest p-value updated : %s" % str(min_p)

  if pval[score+score_range/2]<signif_cutoff:
    if (c1/c2)>(c3/c4): # ((host==1 && word==1)/host==1) > ((host!=1 && word==1)/host!=1)
      whichhost=compare_binary_str_list[0];
    else:
      whichhost=compare_binary_str_list[1];
    lo=-1;
    match=-1;
    for m in range(0,3):
      keys=approxkeys(key,m);
      for key2 in keys:
        if loc.has_key(key2):
          lo=loc[key2];
          match=m;
          break;
      if match>-1:break;
    f1.write(key+','+str(pval[score+score_range/2])+','+str(score)+','+str(scoreOR)+','+whichhost+','+str(int(c1))+','+str(int(c3))+','+str(int(c2-c1))+','+str(int(c4-c3))+','+str(lo)+','+str(match)+','+patstr+'\n');


  if (scoreOR < cutoff_small_scoreOR) | (scoreOR > cutoff_large_scoreOR):
  #if (scoreOR < numpy.min(scoresOR)) | (scoreOR > numpy.max(scoresOR)):

    if (c1/c2)>(c3/c4): # ((host==1 && word==1)/host==1) > ((host!=1 && word==1)/host!=1)
      whichhost=compare_binary_str_list[0];
    else:
      whichhost=compare_binary_str_list[1];
    lo=-1;
    match=-1;
    for m in range(0,3):
      keys=approxkeys(key,m);
      for key2 in keys:
        if loc.has_key(key2):
          lo=loc[key2];
          match=m;
          break;
      if match>-1:break;
    #fOR.write(key+','+str(pval[score+score_range/2])+','+str(score)+','+str(scoreOR)+','+whichhost+','+str(int(c1))+','+str(int(c3))+','+str(int(c2-c1))+','+str(int(c4-c3))+','+str(lo)+','+str(match)+','+patstr+'\n');


f1.close()
#fOR.close()
print 'All done!'


#
# in the future
#   extention to continuous variable
#   how to summarize results of different groups 
#