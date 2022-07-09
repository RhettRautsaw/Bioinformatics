#!/usr/bin/env python

import argparse
import subprocess

parser = argparse.ArgumentParser(description='Takes an RSEM results output and the input fasta and outputs a fasta of the most highly expressed sequences. User can specify how many sequences they want to retain.')
parser.add_argument("-ro","--rsem_output",
					type=str,
                    default='',
                    help="RSEM results sheet based on running RSEM on a file of all ORFs")
parser.add_argument("-i","--input_fasta",
					type=argparse.FileType('r+'),
					help="Fasta of all ORFs")
parser.add_argument("-fdb","--filterdatabase",
					type=argparse.FileType('r+'),
					help="Annotated transcriptome (fasta) used to filter highly expressed ORFs")
parser.add_argument("-db","--blast_database",
					type=str,
                    default="",
                    help="Path to blast database (e.g. SWISSprot)")
parser.add_argument("-remote",
					action="store_true",
                    help="Use remote blast against nucleotide database instead of a local blast database. WARNING: This will take much much much much much longer!")
parser.add_argument("-ol","--output_label",
					type=str,
                    default="",
                    help="Identifier prefix used for all output files")
parser.add_argument("-n","--number",
					type=int,
                    default=500,
                    help="Number of highly expressed ORFs to retain and investigate")
parser.add_argument("-m","--makeblastdb",
					type=str,
                    default="makeblastdb",
                    help="directory with makeblastdb command. Default assumes it is in your path")
parser.add_argument("-bn","--blastn",
					type=str,
                    default="blastn",
                    help="directory with blastn command. Default assumes it is in your path")
parser.add_argument("-p","--matchpercent",
                    nargs='?',
                    type=float,
                    default=99,
                    help="Percent identity for blastn. Default is 99")
parser.add_argument("-bx","--blastx",
					type=str,
                    default="blastx",
                    help="directory with blastx command. Default assumes it is in your path")
parser.add_argument("-nt","--num_threads",
					type=int,
                    default=8,
                    help="Number of threads for blast")
parser.add_argument("-rsc","--RSEM_results_collector",
					type=str,
                    default="~/Dropbox/bin/AndrewScripts/Rogue_CDS_Finder/Step1_RSEM_results_collector.py",
                    help="directory for the 'Step1_RSEM_results_collector.py' script. Default assumes it is in your path")
parser.add_argument("-sif","--Seq_IDnFilter",
					type=str,
                    default="~/Dropbox/bin/AndrewScripts/Rogue_CDS_Finder/Step2_Seq_IDnFilter_blast.py",
                    help="directory for the 'Step2_Seq_IDnFilter_blast.py' script. Default assumes it is in your path")
parser.add_argument("-bnf","--BLAST_nFilter",
					type=str,
                    default="~/Dropbox/bin/AndrewScripts/Rogue_CDS_Finder/Step3_BLAST_nFilter.py",
                    help="directory for the 'Step3_BLAST_nFilter.py' script. Default assumes it is in your path")
args=parser.parse_args()



## Define variables to make my life easier
##########################################################################################
rsc = args.RSEM_results_collector
rsem_output = args.rsem_output
input_fasta = args.input_fasta
input_fasta_name = input_fasta.name
number = args.number
if len(args.output_label) == 0:
	label = args.output_label
else:
	label = args.output_label + '_'

sif = args.Seq_IDnFilter
filterdatabase = args.filterdatabase.name
makeblastdb = args.makeblastdb
blastn = args.blastn
matchpercent = args.matchpercent
num_threads = args.num_threads

bnf = args.BLAST_nFilter
blast_database = args.blast_database
blastx = args.blastx
##########################################################################################

##########################################################################################
print('Collecting top ' + str(number) + ' most highly expressed ORFs')
command = 'python ' + rsc + ' -ro ' + rsem_output + ' -f ' + input_fasta_name + ' -n ' + str(number) + ' -o ' + label + '_hiORF.fasta'
subprocess.call(command, shell=True)
##########################################################################################

##########################################################################################
print('Filtering ORFs using blastn against ' + str(filterdatabase) + ' with a percent identity of ' + str(matchpercent))
command = 'python ' + sif + ' -f ' + label +'_hiORF.fasta -fdb ' + filterdatabase + ' -m ' + makeblastdb + ' -bn ' + blastn + ' -p ' + str(matchpercent) + ' -nt ' + str(num_threads) + ' -o ' + label +'_hiORF_passfilter.fasta'
subprocess.call(command, shell=True)
##########################################################################################

##########################################################################################
print('Blasting ORFs which passed filter and writing output')

if args.remote:
	command = 'python ' + bnf + ' -f ' + label +'_hiORF_passfilter.fasta -remote -x ' + label + '_hiORF_blast.xml -nt ' + str(num_threads) + ' -o ' + label
	subprocess.call(command, shell=True)


else:
	command = 'python ' + bnf + ' -f ' + label +'_hiORF_passfilter.fasta -db ' + blast_database + ' -bx ' + blastx + ' -x ' + label + '_hiORF_blast.xml -nt ' + str(num_threads) + ' -o ' + label
	subprocess.call(command, shell=True)


##########################################################################################

print('Done.')