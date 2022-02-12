#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import sys, os, shutil
import subprocess as sp
import glob
import datetime as dt
from tqdm import tqdm
from Bio import SeqIO
from Bio import AlignIO
from Bio.Nexus import Nexus
from Bio.Alphabet import IUPAC

############################################### ARGUMENTS

parser = argparse.ArgumentParser(description="""
Concatenate a directory of alignments for phylogenetics
""")
parser.add_argument("-i","--input",
					type=str,
					help="directory containing alignments for concatenation")
parser.add_argument("-t","--type",
					type=str,
					default="fasta",
					help="filetype in the directory (e.g. fasta, nexus, phylip). Default is fasta")
parser.add_argument("-o","--output",
					type=str,
					default="concatenate",
                    help="Output directory name")
parser.add_argument("-d","--delimiter",
					type=str,
					default=".",
					help="The delimiter that separates the unique part of a file from the non unique parts. For instance if your files are labelled by gene (e.g. CytB.fas, ND4.fas) the delimiter is . whereas if you files are L1_final.fas, L2_final.fas the delimiter is the _ . A caveat is that this script always assumes the unique part of the name is the first part of the name.")
args=parser.parse_args()

############################################### SETUP

input=os.path.abspath(args.input)
type=args.type
output=args.output
delim=args.delimiter

############################################### CODE

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting alignments to NEXUS and creating partitions :::")
files=sorted(glob.glob(input+'/*'))
sp.call("mkdir -p " + output + "/nexus",shell=True)
x=1
partitions = open(output + "/partitions.txt","w")
for file in tqdm(files):
	filename = os.path.basename(file)
	lociname = filename.split(delim)[0]
	outfile = output + '/nexus/' + lociname + '.nex'
	AlignIO.convert(file, type, outfile,'nexus',alphabet=IUPAC.ambiguous_dna)
	aln = AlignIO.read(file, type)
	partitions.write(lociname + " = " + str(x) + "-" + str(x - 1 + aln.get_alignment_length())+";")
	partitions.write("\n")
	x = x + aln.get_alignment_length()

partitions.close()

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Concatenating NEXUS alignments :::")
loci=sorted([os.path.basename(x).split(".nex")[0] for x in glob.glob(output + "/nexus/*.nex")])
nexi =  [(locus, Nexus.Nexus(output + "/nexus/" + locus + ".nex")) for locus in loci]
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open(output + "/concat.nex", 'w'))


