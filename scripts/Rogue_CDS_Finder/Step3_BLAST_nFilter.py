#!/usr/bin/env python

import argparse
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(description='takes an input fasta and blasts it against a specified database')
parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="Fasta to rename sequences in.")
parser.add_argument("-bx","--blastx",
					type=str,
                    default="blastx",
                    help="directory with blastx command. Default assumes it is in your path")
parser.add_argument("-db","--blast_database",
					type=str,
                    default="",
                    help="path to blast database")
parser.add_argument("-x","--blast_xml",
					type=str,
                    default='',
                    help="name of blast xml file")
parser.add_argument("-remote",
					action="store_true",
                    help="Use remote blast against nucleotide database instead of a local blast database. WARNING: This will take much much much much much longer!")
parser.add_argument("-nt","--num_threads",
					type=int,
                    default=8,
                    help="number of threads for blast")
parser.add_argument("-o","--output",
					type=str,
                    default='',
                    help="conserved part of name for output files (a 'no_hits.fasta' and a 'hits.fasta and hits.csv')")
args=parser.parse_args()


fasta = args.fasta
fasta_name = args.fasta.name
blastx = args.blastx
blast_database = args.blast_database
blast_xml = args.blast_xml
num_threads = args.num_threads
output = args.output
if len(output) == 0 : 
	output = output
else:
	output = args.output + '_'


#fasta = "Remaining_Baurifer_CDS.fasta"
#blast = 'blastx'
#blast_database = 'SWISS-PROT/SWISS-PROT'
#blast_xml = 'Remaing_CDS_blast.xml'
#num_threads = 4


sequences = list(SeqIO.parse(fasta,"fasta"))

if args.remote:
	command = 'blastn -query ' + fasta_name + ' -db nt -remote -outfmt 5 -max_target_seqs 10 -evalue 0.0001 -out ' + blast_xml
	subprocess.call(command, shell=True)

else:
	command = blastx + ' -query ' + fasta_name + ' -db ' + blast_database + ' -outfmt 5 -num_threads ' + str(num_threads) + ' -max_target_seqs 10 -evalue 0.0001 -out ' + blast_xml
	subprocess.call(command, shell=True)

blast_results = SearchIO.parse(blast_xml, 'blast-xml')

#blast = []
#for result in blast_results:
	#blast.append(result)


blast_w_hits = []
blast_wo_hits = []
for rec in blast_results:
	if len(rec) == 0:
		blast_wo_hits.append(rec)
	else:
		blast_w_hits.append(rec)
		

seqs_w_hits = []
for seq in sequences:
	for rec in blast_w_hits:
		if seq.name == rec.id:
			seqs_w_hits.append(seq)


seqs_wo_hits = []
for seq in sequences:
	for rec in blast_wo_hits:
		if seq.name == rec.id:
			seqs_wo_hits.append(seq)


handle=open(output + 'no_hits.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(seqs_wo_hits)
handle.close()

handle=open(output + 'hits.fasta', "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(seqs_w_hits)
handle.close()		

csv_output = []
for rec in blast_w_hits:
	line =[]
	line.append(rec.id)
	line.append(rec.seq_len)
	for hit in rec:
		line.append(hit.description)
	csv_output.append(line)

for rec in blast_wo_hits:
	line = []
	line.append(rec.id)
	line.append(rec.seq_len)
	csv_output.append(line)

with open(output + 'hits.csv', 'w') as csv_file:
	csv_writer = csv.writer(csv_file, delimiter = ',')
	for row in csv_output:
		csv_writer.writerow(row)
		

		
csv_file.close()







