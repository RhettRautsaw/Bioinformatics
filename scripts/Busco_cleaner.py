#!/usr/bin/env python

import sys, os, shutil
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

###############################################

parser = argparse.ArgumentParser(description="""
BUSCO duplicate/isoform cleaner (get more single copy genes)
""")

###############################################

parser.add_argument("-m","--multi",
					type=str,
					default="multi_copy_busco_sequences",
					help="multi_copy_busco_sequences folder")
parser.add_argument("-s","--single",
					type=str,
					default="single_copy_busco_sequences",
					help="single_copy_busco_sequences folder")
args=parser.parse_args()

###############################################

multi = os.path.abspath(args.multi)
single = os.path.abspath(args.single)

multi = os.path.abspath("multi_copy_busco_sequences")
single = os.path.abspath("single_copy_busco_sequences")

single_count = len(glob.glob(single+"/*.fna"))

os.chdir(multi)
if not os.path.exists('Duplicates'):
	os.mkdir('Duplicates')
multi_list = glob.glob("*.fna")
multi_count = len(multi_list)

for locus in multi_list:
	#locus="101012at32523.fna"
	sequences = list(SeqIO.parse(locus,"fasta"))
	
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
	
	if len(unique_seqs) == 1:
		handle=open(single + "/" + locus, "w")
		writer = FastaIO.FastaWriter(handle, wrap=None)
		writer.write_file(unique_seqs)
		handle.close()
		
		shutil.move(locus, 'Duplicates/')
		locus2 = locus.split(".fna")[0] + ".faa"
		shutil.move(locus2, 'Duplicates/')

single_count2 = len(glob.glob(single+"/*.fna"))
multi_count2 = len(glob.glob(multi+"/*.fna"))
dup_count = len(glob.glob(multi+"/Duplicates/*.fna"))

print("Single Copy Busco Loci: "+ str(single_count))
print("Multi Copy Busco Loci: "+ str(multi_count))
print("::::::::::::::")
print("Duplicates in Multi Copy Busco Loci: "+ str(dup_count))
print("::::::::::::::")
print("New Single Copy Busco Loci: "+ str(single_count2))
print("New Multi Copy Busco Loci: "+ str(multi_count2))

