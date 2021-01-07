#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='SortFasta.py\n\n\
This script will sort your fasta by sequence ids or length\n\n\
::CITE:: \n\
https://github.com/RhettRautsaw/Bioinformatics\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta file with duplicated names")
parser.add_argument("-w","--wrap",
					action="store_true",
					default=False,
					help="Default is to linearize, but this flag will wrap your fastas for you.")
parser.add_argument("-l","--length",
					action="store_true",
					default=False,
					help="Default is to sort by sequence ID, but this flag will sort by length instead.")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name
outfile = fasta_name.split('.')[0]
outfile = outfile + "_sorted.fasta"

###############################################

if args.length==True:
	len_and_ids = sorted(
		(len(rec), rec.id) for rec in SeqIO.parse(fasta_name, "fasta")
	)
	ids = reversed([id for (length, id) in len_and_ids])
	del len_and_ids  # free this memory
if args.length==False:
	ids = sorted(
		(rec.id) for rec in SeqIO.parse(fasta_name, "fasta")
	)

record_index = SeqIO.index(fasta_name, "fasta")
records = (record_index[id] for id in ids)

if args.wrap==False:
	handle=open(outfile, "w")
	writer = FastaIO.FastaWriter(handle, wrap=None)
	writer.write_file(records)
	handle.close()

if args.wrap==True:
	handle=open(fasta_name, "w")
	writer = FastaIO.FastaWriter(handle)
	writer.write_file(fasta)
	handle.close()