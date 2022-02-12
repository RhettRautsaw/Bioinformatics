#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""

RemoveSeqs.py

This script will cycle through a fasta file and search for a given list of sequence IDs anc remove them from the fasta.

::CITE:: \n\
https://github.com/RhettRautsaw/Bioinformatics
""")

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta file")
parser.add_argument("-l","--list",
					type=argparse.FileType('r+'),
					help="List (1 name per line) of sequences to remove")
parser.add_argument("-o","--out",
					type=str,
					default=None,
					help="Name of new fasta file")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name
list_name = args.list.name
if args.out is None:
	out_name = fasta_name.split(".fasta")[0]+"_RmSeq.fasta"
else:
	out_name = args.out

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))
list=open(list_name,"r").read().split('\n')

new_fasta=[]
for seq in fasta:
	if seq.name not in list:
		new_fasta.append(seq)

handle=open(out_name, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(new_fasta)
handle.close()

