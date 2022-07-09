#!/usr/bin/env python

import argparse
import csv
import subprocess
from Bio import SeqIO
from Bio.SeqIO import FastaIO

# Command line options
parser = argparse.ArgumentParser(description='This script takes a fasta input and uses cd-hit-est to filter out those represented in the database above a specified sequence identity. For example, it was designed to take a list of CDS\'s and find how many are already represented in an annotated transcriptome (here the filter database)')
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta of CDS's.")
parser.add_argument("-fdb","--filterdatabase",
					type=argparse.FileType('r+'),
					help="faasta database used to filter the highly expressed transcripts. For instance, of your most highly expressed CDS's, how many are already present in a transcriptome")					
parser.add_argument("-o","--output",
					type=str,
                    default=" ",
                    help="name of the output fasta file")
parser.add_argument("-c","--cdhitest",
                    nargs='?',
                    type=str,
                    default="cd-hit-est",
                    help="Path to cd-hit-est executable. Default assumes it is in your PATH.")
parser.add_argument("-p","--matchpercent",
                    nargs='?',
                    type=float,
                    default=0.98,
                    help="Match percentage for cd-hit-est. Default is 0.98")          
args=parser.parse_args()



fasta = args.fasta
filterdatabase = args.filterdatabase
output = args.output
cdhitest = args.cdhitest
match = args.matchpercent


#fasta = 'Bauri-CLP2481_hiCDS.fasta'
#filterdatabase = 'Baurifer.fasta'
#output = 'Remaining_Baurifer_CDS.fasta'
#cdhitest = 'cd-hit-est'
#match = float(0.98)

# Read in Fasta file and database file
sequences = list(SeqIO.parse(fasta,"fasta"))
database = list(SeqIO.parse(filterdatabase,"fasta"))

for seq in database:
	seq.name = "DB__" + seq.name
	seq.id = "DB__" + seq.id
	seq.description = "DB__" + seq.description

cdhit_tmp = []
for seq in sequences:
	cdhit_tmp.append(seq)
	

for rec in database:
	cdhit_tmp.append(rec)
	
	
handle=open("tmp.fasta", "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(cdhit_tmp)
handle.close()

#Run cd-hit to match sequences. 
command = cdhitest + " -i tmp.fasta  -o Clustered.fasta -c " + str(match) + " -T 8 -n 10 -d 0 -M 0 -p 1"
subprocess.call(command,shell=True)

clustered = open('Clustered.fasta.clstr',"r")
clustered = clustered.read()
clustered = clustered.split('>Cluster ')[1:]

clustered_array = []
for line in clustered:
	clustered_array.append(line.split('\n')[:-1])


clustered_CDS = []

for cluster in clustered_array:
	check = 0
	if len(cluster) > 2:
		for seq in cluster:
			if "DB__" in seq:
				check += 1
		for i in range(1,len(cluster)):
			if check > 0 :
				name = cluster[i].split(', >')[1].split('...')[0]
				for seq in sequences:
					if name == seq.name :
						clustered_CDS.append(seq.name)
			if check == 0 :
				if '*' not in cluster[i]:
					name = cluster[i].split(', >')[1].split('...')[0]
				for seq in sequences:
					if name == seq.name :
						clustered_CDS.append(seq.name)


remaining_CDS = [x for x in sequences if x.name not in clustered_CDS]

 
handle=open(output, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(remaining_CDS)
handle.close()			

subprocess.call('rm tmp.fasta',shell=True)
subprocess.call('rm Clustered.fasta.clstr',shell=True)
subprocess.call('rm Clustered.fasta',shell=True)


