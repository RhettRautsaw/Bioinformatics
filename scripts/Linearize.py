#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
Linearize your damn fastas...who uses multi-line fastas...you people are monsters. 

...oh wait...

Yep, multi-line fastas are better when working with giant genomes...alright... -w will wrap those up for you.

::CITE:: 
https://github.com/RhettRautsaw/Bioinformatics\n\n
""")

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="fasta")
parser.add_argument("-w","--wrap",
					action="store_true",
					default=False,
					help="Default is to linearize, but this flag will wrap your fastas for you.")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))

if args.wrap==False:
	handle=open(fasta_name, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(fasta)
	handle.close()

if args.wrap==True:
	handle=open(fasta_name, "w")
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(fasta)
	handle.close()