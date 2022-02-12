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
					type=str,
					help="Fasta to rename sequences in.")
parser.add_argument("-w","--wrap",
					action="store_true",
					default=False,
					help="Default is to produce linear fasta, but this flag will wrap your fastas for you.")
parser.add_argument("-o","--output",
					type=str,
					default="",
					help="Fasta to rename sequences in.")
args=parser.parse_args()

if args.output=="":
	output=args.fasta
else:
	output=args.output

## Parse sequences, store as list
sequences = list(SeqIO.parse(args.fasta,"fasta"))

for i in sequences:
	i.seq = i.seq.ungap("-")

## Write results to file
if args.wrap==False:
	SeqIO.write(sequences, output, "fasta-2line")
if args.wrap==True:
	SeqIO.write(sequences, output, "fasta")
