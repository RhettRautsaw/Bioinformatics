#!/usr/bin/env python

# Requirements:
#	- Biopython
#	- PyVCF

import sys, os, shutil
import datetime as dt
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess as sp
import gzip
import re
from tqdm import tqdm
try:
	import vcf
except:
	print("Error: PyVCF module is not properly installed.")
	quit()

try:
	from Bio import SeqIO
	from Bio.SeqIO import FastaIO
	from Bio import AlignIO
	from Bio.Nexus import Nexus
except:
	print("Error: biopython module is not properly installed.")
	quit()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
Script will read a VCF and convert it to a concatenated fasta and nexus for input into SVDQuartets. 
It will produce separate sequences for each allele (i.e., Sample1_a1, Sample1_a2)
and will include a partition to coalesce the two alleles together in SVDQuartets

:: EXAMPLE ::
ConcatSNPs.py -i example.vcf

:: CITE :: 
https://github.com/RhettRautsaw/Bioinformatics\n\n""")

###############################################

parser.add_argument("-i","--input",
					type=str,
					default=None,
					help="VCF file containing SNPs")
args=parser.parse_args()

############################################### SETUP
#input=os.path.abspath("tumors.recode.vcf")
input=os.path.abspath(args.input)
name=input.split("/")[-1].split(".vcf")[0]

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Reading VCF :::")
vcf_reader=vcf.Reader(open(input,'r'))

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Creating Concatenated Dictionary :::")
sample_dict=dict()
sample_dict2=dict()
for record in tqdm(vcf_reader):
	for sample in record.samples:
		key0=sample.sample.replace("-","_")
		key1=sample.sample.replace("-","_")+"_1"
		key2=sample.sample.replace("-","_")+"_2"
		if key1 not in sample_dict:
			sample_dict[key1]=list()
			sample_dict[key2]=list()
		if sample.gt_bases != None:
			sample_dict[key1].extend(re.split('/|\|', sample.gt_bases)[0])
			sample_dict[key2].extend(re.split('/|\|', sample.gt_bases)[1])
		else:
			sample_dict[key1].extend('-')
			sample_dict[key2].extend('-')
		if key0 not in sample_dict2:
			sample_dict2[key0]=[key1,key2]

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting to Fasta :::")
ofile = open(name + ".fasta", "w")
for sample in sample_dict:
	concat=''.join(sample_dict[sample])
	ofile.write(">" + sample + "\n" + concat + "\n")

ofile.close()

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting to Nexus :::")
AlignIO.convert(name+".fasta", "fasta", name+".nex", 'nexus', "DNA")
n=Nexus.Nexus(name+".nex")
n.taxpartitions={"samples": sample_dict2}
n.taxsets={"samples": sample_dict2}
n.write_nexus_data(name+".nex", interleave=False, append_sets=True, codons_block=False)

sp.call("perl -pi -e 's/^taxset.*\n//g' " + name + ".nex", shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Finished :::")