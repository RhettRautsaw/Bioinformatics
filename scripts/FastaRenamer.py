#!/usr/bin/env python

import argparse
import random
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO


# Command line options
parser = argparse.ArgumentParser(description='')
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta to rename sequences in.")
parser.add_argument("-n","--name",
					type=str,
					default="",
					help="prefix you want added to all sequences")
parser.add_argument("-w","--wrap",
					action="store_true",
					default=False,
					help="Default is to produce linear fasta, but this flag will wrap your fastas for you.")
args=parser.parse_args()


## Parse sequences, store as list
sequences = list(SeqIO.parse(args.fasta,"fasta"))

for i in range (0,len(sequences)):
	sequences[i].id = args.name + str(i)
	sequences[i].description = ""
		
args.fasta.seek(0)
args.fasta.truncate()
## Write results to file
## Line below hashed out in V2 to write sequences as a single line
#SeqIO.write(best_seqs, args.fasta, "fasta")
if args.wrap==False:
	fasta_out = FastaIO.FastaWriter(args.fasta, wrap=None)
	fasta_out.write_file(sequences)
if args.wrap==True:
	fasta_out = FastaIO.FastaWriter(args.fasta)
	fasta_out.write_file(sequences)
