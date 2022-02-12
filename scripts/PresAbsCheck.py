#!/usr/bin/env python

# Additional software necessary to run this:
# (1) bwa 
# (2) samtools
# (3) bedtools
# (4) biopython
# (5) pandas
# (6) dfply

import argparse
import sys, os, shutil
import subprocess as sp
import datetime as dt
import csv
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from dfply import *

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(description='Check for transcripts with low coverage across a certain percentage of the contig/transcript. This may be indicitive of a misassembly or false positive annotation. Default is to check for transcripts with <5x coverage for >10% of the total contig/transcript length')
parser.add_argument("-i","--input",
					type=str,
					help="Fasta file of contigs/transcripts to check. Use only CODING SEQUENCES.")
parser.add_argument("-r","--reads",
					type=str,
					help="Fastq file of UNPAIRED reads.")
parser.add_argument("-c","--cov",
					type=int,
					default=5,
					help="Coverage threshold. (default: %(default)s)")
parser.add_argument("-p","--prop",
					type=int,
					default=10,
					help="Proportion of the contig/transcript. (default: %(default)s)")
parser.add_argument("-t","--threads",
					type=int,
					default=1,
					help="Number of processing threads. (default: %(default)s)")
args = parser.parse_args()

########################################
################# SETUP ################
########################################

input = os.path.abspath(args.input)
input_name = os.path.basename(input)
input_name2 = input_name.split(".fasta")[0]
reads = os.path.abspath(args.reads)
cov = args.cov
prop = args.prop
threads = args.threads

print("""

Transcript Presence/Absence Checker

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Starting...")
print("\tInput --> "+ input_name)
print("\tReads --> "+ os.path.basename(reads))
print("\tRemoving transcripts with...")
print("\t\t< " + str(cov) + "x coverage")
print("\t\tfor > " + str(prop) + "% of the total transcript length")
print("\tThreads --> " + str(threads))

########################################
################# CODE #################
########################################

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Generating bwa alignment :::")
sp.call("bwa index " + input, shell=True)
sp.call("bwa mem -t " + str(threads) + " " + input + " " + reads + " | samtools sort -@ " + str(threads) + " > " + input_name2 + ".bam", shell=True)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Calculating coverage :::")
sp.call("bedtools genomecov -d -ibam " + input_name2 + ".bam > coverage.txt", shell=True)
sp.call("rm " + input + ".* " + input_name2 + ".bam", shell=True)

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Assessing coverage for Presence/Absence :::")
coverage = pd.read_csv("coverage.txt",sep='\t',names=['transcript','pos','cov'])
seqs = list(SeqIO.parse(input,'fasta'))
print("Input transcripts = " + str(len(seqs)))
print("Coverage data for = " + str(len(set(coverage['transcript']))))

tmp1 = coverage >>  group_by(X.transcript) >> summarize(length=n(X.pos))
tmp2 = coverage[coverage['cov'] < cov]
if(len(tmp2)!=0):
	tmp2 = tmp2 >> group_by(X.transcript) >> summarize(lowcov_count=n(X.pos))
	tmp3 = tmp2 >> right_join(tmp1, by='transcript')
	tmp4 = tmp3 >> mutate(lowcov_prop = (X.lowcov_count/X.length)*100)
	tmp4.to_csv("coverage_low.txt", sep="\t", index=False)
	tmp5 = tmp4[tmp4['lowcov_prop'] > prop]
	T = tmp5['transcript'].tolist()
else:
	T = list()

print(str(len(T)) + " transcripts with less than " + str(cov) + "x coverage for greater than " + str(prop) + "% of total transcript length")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Removing \"Absent\" sequences :::")

good_seqs = []
for seq in seqs:
	if seq.name not in T:
		good_seqs.append(seq)

output_handle=open(input_name2 + "_present.fasta", "w")
writer = FastaIO.FastaWriter(output_handle, wrap=None)
writer.write_file(good_seqs)
output_handle.close()