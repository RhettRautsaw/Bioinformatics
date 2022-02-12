#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import datetime as dt
import random
import sys, os, shutil
import subprocess as sp
import glob

# Command line options
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
fastq-dump and fasterq-dump are dumb and have dumb defaults. Let me help format those sequences and change the quality scores for you into something that makes sense and won't screw with downstream data processing. 

Just give me: 
1. The SRA/SRX/SRR number you want
2. A prefix for renaming files

I will grab the data, name the reads so that it doesn't screw with downstream processing, gzip as quickly as possible, concatenate SRR runs together, and reformat quality scores to PHRED+33 (i.e., the normal one).

I'm still pretty slow, but it's not my fault. It's fastq-dump's fault...Did I mention that fastq-dump is dumb?
""")
parser.add_argument("-s","--sra",
					type=str,
					default=None,
					help="SRA/SRX/SRR number to pull")
parser.add_argument("-n","--name",
					type=str,
					default=None,
					help="Prefix you want added to all files")
parser.add_argument("--fasterq",
					action="store_true",
					default=False,
					help="Use fasterq-dump instead of fastq-dump?")
parser.add_argument("-t" , "--threads",
					type=int,
					default=8,
					help="I'll try to parallelize as much as I can. How many threads should I use?")
args=parser.parse_args()

if args.sra==None:
	print("Come on...you gotta at least give me a SRA number to work with....get your shit together man...")
	quit()

os.mkdir(args.sra)
os.chdir(args.sra)

if args.fasterq:
	sp.call("fasterq-dump -p -x -e " + str(args.threads) + " " + args.sra, shell=True)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Renaming fastq Sequences :::")
	files=glob.glob("*.fastq")
	for i in files:
		j=i.replace("_",".").split(".")[0]
		sp.call('seqtk rename ' + i + ' "' + j + ':" > tmp.fq', shell=True)
		sp.call('mv tmp.fq ' + i, shell=True)
	#sp.call("perl -pi -e 's/\+.*/+/g' *.fastq", shell=True)
	#sp.call("perl -pi -e 's/^(@.*)\.(.*?) .*/$1:$2/g' *.fastq", shell=True)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Adding Read IDs :::")
	sp.call("perl -pi -e 's/^(@.*?) .*/$1 1/g' *_1.fastq", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	sp.call("perl -pi -e 's/^(@.*?) .*/$1 2/g' *_2.fastq", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Gzipping fastq Sequences :::")
	sp.call("pigz -9 -p " + str(args.threads) + " *.fastq", shell=True)

	if args.name!=None:
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Concatenating SRR Runs :::")
		sp.call("cat *_1.fastq.gz > " + args.name + "_R1.fastq.gz", shell=True,stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call("cat *_2.fastq.gz > " + args.name + "_R2.fastq.gz", shell=True,stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call('cat "$(ls *.fastq.gz | grep -v _)" > ' + args.name + '_R0.fastq.gz', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call('rm *_1.fastq.gz *_2.fastq.gz "$(ls *.fastq.gz | grep -v _)"', shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

else:
	sp.call("prefetch " + args.sra, shell=True)
	sp.call("fastq-dump --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --defline-seq '@$ac:$si $ri' --defline-qual '+' */*.sra", shell=True)
	print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Gzipping fastq Sequences :::")
	sp.call("pigz -9 -p " + str(args.threads) + " *.fastq", shell=True)
	sp.call("rm -r */", shell=True)

	if args.name!=None:
		print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Concatenating SRR Runs :::")
		sp.call("cat *_pass_1.fastq.gz > " + args.name + "_R1.fastq.gz", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call("cat *_pass_2.fastq.gz > " + args.name + "_R2.fastq.gz", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call("cat *_pass.fastq.gz > " + args.name + ".fastq.gz", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		sp.call("rm *_pass*.fastq.gz", shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Removing Empty Files :::")
sp.call("find . -size  0 -print -delete", shell=True)

print(dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" ::: Reformatting Quality Scores :::")
files=glob.glob("*.fastq.gz")
for i in files:
	j=i.split(".fastq.gz")[0]
	sp.call("reformat.sh in=" + i + " out=" + j + ".2.fastq.gz qout=33 overwrite=true", shell=True)
	sp.call("mv " + j + ".2.fastq.gz " + i, shell=True)

sp.call("mv *.fastq.gz ../", shell=True)
os.chdir("..")
os.rmdir(args.sra)