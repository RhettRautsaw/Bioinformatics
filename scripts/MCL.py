#!/usr/bin/env python3

import sys, os, shutil, errno
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
#import multiprocessing as mp
#import glob
#from dfply import *
#from p_tqdm import p_map
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqIO import FastaIO
except:
	print("Error: biopython module is not properly installed.")
	quit()

########################################
############### ARGUMENTS ##############
########################################

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

Cluster sequences and keep sequences in largest cluster as reference

""")

parser.add_argument("-i","--input",
					type=str,
					help="Input fasta")
parser.add_argument("--inflation",
					type=float,
					default=2.0,
					help="MCL Inflation Parameter (default: %(default)s)")
parser.add_argument("-c","--cpu",
					type=int,
					default=mp.cpu_count(),
					help="Number of threads to be used in each step. (default: %(default)s)")
args=parser.parse_args()

########################################
############### FUNCTIONS ##############
########################################

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

########################################
################# SETUP ################
########################################

input = os.path.abspath(args.input)
input_name = os.path.basename(input)
input_name2 = input_name.split(".fasta")[0]

inflation = args.inflation
cpus = args.cpu

print("""

Cluster sequences

""")

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: starting BLAST and MCL...")
print("\tInput Fasta -> "+ input)
print("\tMCL Inflation Parameter -> " + str(inflation))
print("\tNumber of CPU -> "+ str(cpus))

########################################
################# CODE #################
########################################

mkdir_p("MCL_results")

seqs = list(SeqIO.parse(input,'fasta'))

print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running BLAST :::")
sp.call("makeblastdb -in " + input + " -dbtype nucl -out MCL_results/" + input_name, shell=True)
sp.call("blastn -query " + input + " -db MCL_results/" + input_name + " -evalue 1e-10 -task dc-megablast -outfmt 6 -out MCL_results/" + input_name + ".blast", shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Running MCL :::")
os.chdir("MCL_results")
sp.call("cut -f 1,2,11 " + input_name + ".blast > " + input_name + ".abc", shell=True)
sp.call("mcxload -abc " + input_name + ".abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + input_name + ".mci -write-tab " + input_name + ".tab", shell=True)
sp.call("mcl " + input_name + ".mci -I " + str(inflation) + " -use-tab " + input_name + ".tab", shell=True)
sp.call("head -1 out." + input_name + ".mci.I* | sed 's/\\t/\\n/g' > " + input_name + ".mcl", shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Subsetting Cluster 1 Sequences :::")
seqs2=[]
clust1 = open(input_name + ".mcl", "r").read().split("\n")
for seq in seqs:
	if seq.name in clust1:
		seqs2.append(seq)

output_handle = open(input_name2 + ".clust.fasta", 'w')
SeqIO.write(seqs2, output_handle, "fasta")
output_handle.close()

sp.call("rm *" + input_name + ".*", shell=True)
