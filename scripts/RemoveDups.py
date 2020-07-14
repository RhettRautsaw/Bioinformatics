#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='RemoveDups.py\n\n\
This script will cycle through a fasta file and search for duplicate sequence IDs. If found, it will keep the first sequence it encounters\n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta file with duplicated names")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name
outfile = fasta_name.split('.')[0]
outfile = outfile + "_noDups.fasta"

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))

new_fasta=[]
names=[]
for seq in fasta:
	if seq.id not in names:
		names.append(seq.id)
		new_fasta.append(seq)

handle=open(outfile, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(new_fasta)
handle.close()
