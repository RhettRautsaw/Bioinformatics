#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='Removes the last three nucleotides (i.e. stop codon) from sequences. \n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="fasta")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name
outfile = fasta_name.split('.')[0]
outfile = outfile + "_noStop.fasta"

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))

for seq in fasta:
	seq.seq = seq.seq[0:-3]

handle=open(outfile, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(fasta)
handle.close()
