#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
#from sets import Set

# Command line options
parser = argparse.ArgumentParser(description='Some dumb f*cked script I am writing')
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta to rename sequences in.")
parser.add_argument("-o","--out",
					type=str,
                    default="",
                    help="will make a _duplicates.fasta and a _ambiguities.fasta based on the name (e.g. Catro-CLPXXX_duplicates.fasta")
args=parser.parse_args()

sequences = list(SeqIO.parse(args.fasta,"fasta"))

seq_list = []
for seq in sequences:
	seq_list.append(seq.id)


duplicates = []
dup_string = ''
for seq in sequences:
	for seq1 in sequences:
		if seq.seq == seq1.seq and seq.id != seq1.id and seq1.id not in dup_string:
			duplicates.append(seq)
			dup_string = dup_string + seq.id
			break
		
			
seq_list = [x for x in seq_list if x not in dup_string]


unique_seqs = []
for seq in sequences:
	for seq1 in seq_list:
		if seq.id == seq1:
			unique_seqs.append(seq)


ambiguities = []
ambiguities_ids = []
alphabet = ['A','a','C','c','T','t','G','g']			
for seq in unique_seqs:
	for site in seq.seq:
		if site not in alphabet:
			ambiguities.append(seq)
			ambiguities_ids.append(seq.id)
			break


seq_list = [x for x in seq_list if x not in ambiguities_ids]


good_seqs = []
for seq in unique_seqs:
	for seq1 in seq_list:
		if seq.id == seq1:
			good_seqs.append(seq)


write_good = args.out + "_clean.fasta"
handle=open(write_good, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(good_seqs)
handle.close()

write_dup = args.out + "_duplicates.fasta"
handle=open(write_dup, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(duplicates)
handle.close()

write_ambig = args.out + "_ambiguities.fasta"
handle=open(write_ambig, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(ambiguities)
handle.close()
