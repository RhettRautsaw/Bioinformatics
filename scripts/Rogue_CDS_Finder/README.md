# Rogue CDS Finder

This script is designed help detect highly expressed ORFs/CDSs that may have been missed during annotation or accidentally removed during another processing step (*i.e.*, Rogue CDSs).

With a fasta of all possible ORFs from an assembly and their relative expression, Rogue_CDS_Finder will:
- **Step 1**: Extract the # most highly expressed ORFs. 
- **Step 2**: Compare the high-expression ORFs to an annotated transcriptome and remove any ORFs with a match. This results in a fasta of high-expression CDSs that are potentially novel, unidentified, or otherwise missed previously. 
- **Step 3**: Perform a BLAST search of these Rogue CDSs to a database (*e.g.*, SWISSprot)

We split this script up into three parts (as above) which can be ran separately, or can be ran used together using the wrapper function.

## Arguments
| Input						| Description											|
|---------------------------|-------------------------------------------------------|
| `-ro` <br> `--rsem_output`	| RSEM results sheet based on running RSEM on a file of all ORFs |
| `-i` <br> `--input_fasta`		| Fasta of all ORFs |
| `-fdb` <br> `--filterdatabase`| Annotated transcriptome (fasta) used to filter highly expressed ORFs |
| `-db`<br> `--blast_database`	| Path to blast database (e.g. SWISSprot) |
| `-remote`					| Use remote blast against nucleotide database instead of a local blast database. WARNING: This will take much much much much much longer! |
| `-ol`<br> `--output_label`	| Identifier prefix used for all output files |
| `-n` <br> `--number`			| Number of highly expressed ORFs to retain and investigate [Default: 500] |
| `-p` <br> `--matchpercent`	| Percent identity for blastn [Default: 99] |
| `-nt`<br> `--num_threads`	| Number of threads for blast [Default: 8] |

