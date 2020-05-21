#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='RenameDups.py\n\n\
This script will cycle through a fasta file and search for duplicate sequence IDs. If found, it will add a unique index to the sequence ID\n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta file with duplicated names")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))

names=[]
for idx, seq in enumerate(fasta):
	if seq.id in names:
		names.append(seq.id)
		seq.id=seq.id+"_"+str(idx)
		seq.name=seq.name+"_"+str(idx)
	else:
		names.append(seq.id)

handle=open('out.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(fasta)
handle.close()

