#!/usr/bin/env python

import argparse
import csv
from Bio import SeqIO
from Bio.SeqIO import FastaIO

# Command line options
parser = argparse.ArgumentParser(description='Takes an RSEM results output and the input fasta and outputs a fasta of the most highly expressed sequences. User can specify how many sequences they want to retain.')
parser.add_argument("-ro","--rsem_output",
					type=str,
                    default=" ",
                    help="RSEM results sheet based on running RSEM on a file of all ORFs")
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta of all ORFs")
parser.add_argument("-n","--number",
					type=int,
                    default=500,
                    help="Number of highly expressed ORFs to retain and blast")
parser.add_argument("-o","--output",
					type=str,
                    default=" ",
                    help="Name of the output fasta file")
args=parser.parse_args()

def ReadRSEMResults(X):
	RSEM = []
	file = open(X, "r")
	file1 = file.readlines()
	for line in file1:
		line=line.split('\n')[0]
		RSEM.append(line.split('\t'))
	return RSEM

def TakeSixthItem(element):
	return float(element[5])


rsem_file = args.rsem_output
fasta = args.fasta
output = args.output
number = args.number

#rsem_file = 'Bauri-CLP2481_ORF_results.isoforms.results'
#number = 500
#fasta = 'Bauri-CLP2481_ORFs_clustered.fasta' 
#filter = 'Baurifer.fasta'

rsem = ReadRSEMResults(rsem_file)
header = rsem[0]
rsem = rsem[1:]
rsem.sort(key=TakeSixthItem, reverse = True)
rsem = rsem[0:number]


CDS = list(SeqIO.parse(fasta,"fasta"))

hi_CDS = []
for rec in rsem:
	for seq in CDS:
		if rec[0] == seq.name:
			hi_CDS.append(seq)
			break


handle=open(output, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(hi_CDS)
handle.close()
