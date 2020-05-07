#!/usr/bin/env python

import sys, os, shutil
import argparse
from argparse import RawTextHelpFormatter
import csv
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='Linearize your damn fastas...who uses multi-line fastas...you people are monsters. \n\n\
::CITE:: \n\
https://github.com/reptilerhett/Bioinformatic-Scripts\n\n')

###############################################

parser.add_argument("-f","--fasta",
					type=argparse.FileType('r+'),
					help="fasta")
args=parser.parse_args()

###############################################

fasta_name = args.fasta.name

###############################################

fasta = list(SeqIO.parse(fasta_name,"fasta"))

handle=open(fasta_name, "w")
writer = FastaIO.FastaWriter(handle, wrap=None)
writer.write_file(fasta)
handle.close()
