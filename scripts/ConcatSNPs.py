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
It is polyploid friendly and will produce separate sequences for each allele (i.e., Sample1_a1, Sample1_a2, Sample1_a3, etc.).
and will include a partition to coalesce the alleles together per sample in SVDQuartets

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
allele_dict=dict()
sample_dict=dict()
for record in tqdm(vcf_reader):
	for sample in record.samples:
		key0=sample.sample.replace("-","_")
		alleles=sample.ploidity
		for i in range(alleles):
			key=sample.sample.replace("-","_")+"_"+str(i+1)
			if key not in allele_dict:
				allele_dict[key]=list()
			if sample.gt_bases != None:
				allele_dict[key].extend(re.split('/|\|', sample.gt_bases)[i].replace("*","-"))
			else:
				allele_dict[key].extend('-')
			if key0 not in sample_dict:
				sample_dict[key0]=[key]
			if key not in sample_dict[key0]:
				sample_dict[key0].append(key)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting to Fasta :::")
ofile = open(name + ".fasta", "w")
for sample in allele_dict:
	concat=''.join(allele_dict[sample])
	ofile.write(">" + sample + "\n" + concat + "\n")

ofile.close()

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Converting to Nexus :::")
AlignIO.convert(name+".fasta", "fasta", name+".nex", 'nexus', "DNA")
n=Nexus.Nexus(name+".nex")
n.taxpartitions={"samples": sample_dict}
n.taxsets={"samples": sample_dict}
n.write_nexus_data(name+".nex", interleave=False, append_sets=True, codons_block=False)

sp.call("perl -pi -e 's/^taxset.*\n//g' " + name + ".nex", shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Finished :::")