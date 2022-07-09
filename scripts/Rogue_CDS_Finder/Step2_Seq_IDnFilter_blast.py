#!/usr/bin/env python

import argparse
import csv
import subprocess
from Bio import SeqIO
from Bio.SeqIO import FastaIO

# Command line options
parser = argparse.ArgumentParser(description='This script takes a fasta input and uses blast to filter out those represented in the database above a specified sequence identity. For example, it was designed to take a list of CDS\'s and find how many are already represented in an annotated transcriptome (here the filter database)')
parser.add_argument("-f","--fasta",
					type=str,
					help="Fasta of all ORFs")
parser.add_argument("-fdb","--filterdatabase",
					type=str,
					help="Annotated transcriptome (fasta) used to filter highly expressed ORFs")					
parser.add_argument("-o","--output",
					type=str,
                    default=" ",
                    help="Name of the output fasta file")
parser.add_argument("-m","--makeblastdb",
                    nargs='?',
                    type=str,
                    default="makeblastdb",
                    help="Path to makeblastdb executable. Default assumes it is in your PATH.")
parser.add_argument("-bn","--blastn",
                    nargs='?',
                    type=str,
                    default="blastn",
                    help="Path to blastn executable. Default assumes it is in your PATH.")
parser.add_argument("-p","--matchpercent",
                    nargs='?',
                    type=float,
                    default=99,
                    help="Percent identity for blastn. Default is 99")
parser.add_argument("-nt","--num_threads",
					type=int,
                    default=8,
                    help="number of threads for blast")
args=parser.parse_args()



fasta = args.fasta
filterdatabase = args.filterdatabase
output = args.output
makeblastdb = args.makeblastdb
blastn = args.blastn
match = args.matchpercent
num_threads = args.num_threads



# Read in Fasta file and database file
sequences = list(SeqIO.parse(fasta,"fasta"))
database = list(SeqIO.parse(filterdatabase,"fasta"))

command = makeblastdb + " -dbtype nucl -in " + filterdatabase + " -out tmp.blastdb"
subprocess.call(command,shell=True)

command =  blastn + " -query " + fasta + " -out match.blast.out -db tmp.blastdb -perc_identity " + str(match) + " -outfmt 6 -num_threads " + str(num_threads)
subprocess.call(command,shell=True)

command = "cut -d$'\t' -f1 match.blast.out | uniq > tmp.txt"
subprocess.call(command,shell=True)


with open('tmp.txt', 'r') as f:
	clustered = f.read().splitlines()

remaining_CDS = [x for x in sequences if x.name not in clustered]

handle=open(output, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(remaining_CDS)
handle.close()	

subprocess.call('rm tmp.txt',shell=True)
subprocess.call('rm match.blast.out',shell=True)
subprocess.call('rm tmp.blastdb*',shell=True)