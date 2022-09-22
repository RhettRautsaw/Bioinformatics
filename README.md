# Bioinformatics
## Rhett M. Rautsaw
This repository contains bioinformatic advice along with many simple, but potentially useful personally-developed bioinformatics scripts, functions, pipelines and resources for improving your skillset and performing really cool analyses.

# Table of Contents
- [General Advice/Training](#general-advicetraining)
- [Personally Developed Analysis Resources](#personally-developed-analysis-resources)
	- [Pipelines](#pipelines)
	- [Workshops/Tutorials](#workshopstutorials)
	- [Scripts](#scripts)
- [Additional Resources](#additional-resources)
	- [Misc Functions](#misc-functions)
	- [Links](#links)

# General Advice/Training
If you are interested in learning bioinformatics, Unix, [Python](https://www.python.org/), or R then I recommend you check out the [Basics_of_Bioinformatics.md](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials/Basics_of_Bioinformatics.md) document in the [tutorials](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials) folder in this repository. In this document I have provided a variety of workshops and advice for how to become more proficient, efficient, and productive! Even if you have experience in bioinformatics, make sure to check out the section on [GNU-Parallel](https://www.gnu.org/software/parallel/sphinx.html#)! I love GNU-Parallel so much that I was even given a [GNU-Parallel t-shirt](https://gnuparallel.threadless.com/) as a gift. Nerdy? Yes, but that should also tell you how amazing this program is.

Other than the [Basics_of_Bioinformatics.md](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorial/Basics_of_Bioinformatics.md) document in the [tutorials](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials) folder, there is one piece of advice that I think is critical for being a good bioinformatician regardless of your level of experience...

<center> <h2> STAY ORGANIZED </h2> </center>

This is the \#1 piece of advice I give to everyone in bioinformatics. In fact, I would argue that organization is one of the most important skills you can have as a bioinformatician. So how do you do it? Well there are two main things I can tell you, but you have to find a system that works for you.

First, I **HIGHLY** recommend that you keep well documented READMEs for each step you take and redos of a given analysis. I typically create my READMEs for each project I am working on. I create the READMEs in [RMarkdown](https://rmarkdown.rstudio.com/lesson-1.html) format in [RStudio](https://www.rstudio.com/products/rstudio/). This allows me to create code "chunks" to document `unix`, `python`, and `R` code. I can even run these code chunks from RStudio! I then knit them to a GitHub MarkDown document like this one you're currently reading in order to document my code openly for publication. You can see an example of this [here]().

Second, keep your directories ordered and well-labeled. My directories follow a strict format to monitor what samples have what sequencing data associated with them as well as what has been done for each. I typically number my folders to keep track of the order in which I did things. This pattern carries across all my sequence processing regardless of whether it is RNA-Seq, WGS, or something else.

| Folder | Description | 
|--------|-------------|
| `00_raw` | Raw data + [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [NanoPlot](https://github.com/wdecoster/NanoPlot) quality check and statistics |
| `01_concat` | Concatenated forward and reverse reads or multiple sequencing runs + [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [NanoPlot](https://github.com/wdecoster/NanoPlot) quality check and statistics |
| `02_trim` | Quality trimming with [Trim-Galore](https://github.com/FelixKrueger/TrimGalore) or [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) |
| `03_merge` | Merging reads with [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) |
| `04_[assembler]` | All assemblies are started with the `04_` prefix and placed in their own folder (*e.g.*, [trinity](https://github.com/trinityrnaseq/trinityrnaseq), [extender](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-312), [ngen](https://www.dnastar.com/manuals/seqman-ngen/17.0/en/topic/welcome-to-seqman-ngen), [masurca](https://github.com/alekseyzimin/masurca), [spades](https://github.com/ablab/spades), [hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html), etc.). | 
| `05_final-assembly` | Assemblies with renamed headers (`Fasta_Renamer.py`)) and concatenated assemblies | 
| `06_busco` | [BUSCO](https://busco.ezlab.org/) runs on the final assemblies (typically only run on the Trinity assembly for RNA-seq) |
| `07_[annotation]` | All annotation steps are started with the `07_` prefix and placed in their own folders (*e.g.*, [toxcodan](https://github.com/pedronachtigall/ToxCodAn), [nontoxins](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide), [funannotate](https://funannotate.readthedocs.io/en/latest/index.html), [braker](https://github.com/Gaius-Augustus/BRAKER), etc.). | 
| `08_[cleaning]` | Any cleaning steps done to the final annotated fasta start with the `08_` prefix and are placed in their own folders (*e.g.*, `CDS_filter.py`, [`ChimeraKiller`](https://github.com/masonaj157/ChimeraKiller)) | 
| `09_[mapping]` | Any mapping steps done to a reference for the purpose of calling variants, estimating expression, etc. start with `09_` prefix and are placed in their own folders (*.e.g.*, [bwa](http://bio-bwa.sourceforge.net/), [rsem](https://github.com/deweylab/RSEM), [kallisto](https://pachterlab.github.io/kallisto/about), [minimap](https://github.com/lh3/minimap2), [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), etc.) | 
| `10_[analysis]` | All other analyses such as variant calling, scaffolding, etc. start with `10_` prefix and are place din their own folders |
| `99_mitosis` | Extraction of the mitochondrial genome using [MitoSIS](https://github.com/RhettRautsaw/MitoSIS) to check for contamination |
| `99_finalome` | Folder containing the final assemblies/annotation files |

**Example Directory Structure:**
```
Species-SampleNumber/
└── DataType
	└── TissueType
		└── Date_Location_Sequencer_RunName_ReadLength
			├── 00_raw
			├── 01_concat
			├── 02_trim
			...

Ctigr-CLP2741/
├── RNA
│   └── VG-B
│  		└── 2019-03-15_CU_NextSeq_CLP03_150PE
│  			├── 00_raw
│  			├── 01_concat
│  			├── 02_trim
│  			├── 03_merge
│  			...
└── WGS
	└── blood
		├── 2018-11-15_UDel_PacBio_CLR
		│   ├── 00_raw
		│   ├── 01_concat
		│   ├── 04_masurca
		│	...
		└── 2019-01-08_FSU_NovaSeq_DRR19-DRR20_150PE
			├── 00_raw
			├── 01_concat
			├── 02_trim
			...
```

# Personally-Developed Analysis Resources
## Pipelines
These are pipelines I have developed or co-developed and formally scripted. These can generally be found on my other GitHub repositories.

- [MitoSIS](https://github.com/RhettRautsaw/MitoSIS)
	- Python script used to extract the mitochondrial genome using [MITGARD](https://github.com/pedronachtigall/MITGARD) and compare it to Genbank sequences to look for potential contamination.
- [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn)
	- Python script used to annotate snake venom toxin genes using generalized Hidden Markov Models
	- [Hitchhiker's Guide to Venom Gland Transcriptomics](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide)
		- Generalized pipeline for processing snake venom gland transcriptomic data including trimming, merging, assembly, annotation, chimera removal, clustering, expression estimation, visualization, and analysis. 
- [autokuenm](https://github.com/RhettRautsaw/VenomMaps/tree/master/code/autokuenm)
	- Unix-executable R script used to automatically construct Species Distribution Models (SDMs) with [kuenm](https://github.com/marlonecobos/kuenm) after preparing M-areas and partitioning occurrence records.
- [PhyProbe](https://github.com/RhettRautsaw/PhyProbe)
	- Python script used to extract phylogenomic loci from a variety of different sources such as RNA-Seq, WGS, Sequence Capture, etc. 

## Workshops/Tutorials
These are workshops/tutorials/pipelines I have put together, but have not formally scripted (yet). These may be found here or on my other GitHub repositories depending on the size of the required files and/or intensity of the subject matter. 

- [GIS Tutorial](https://github.com/RhettRautsaw/GIS_Tutorial)
	- [QGIS Tutorial](https://github.com/RhettRautsaw/GIS_Tutorial/blob/master/QGIS_Tutorial.md)
	- [Spatial Data in R](https://github.com/RhettRautsaw/GIS_Tutorial/blob/master/R_Tutorial.md)
- [Basics of Bioinformatics](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials/Basics_of_Bioinformatics.md)
	- The basics of how to use Unix
- [HPC Cliffnotes](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials/HPC_Cliffnotes.md)
	- The basics of how to efficiently use HPCs (PBS & translation to SLURM)
- [PacBio HiFi Genomics](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials/HiFi_Genomics.md)
	- General pipeline for assembling and annotating PacBio HiFi genomic data.
- [TreePL Divergence Dating](https://github.com/RhettRautsaw/Bioinformatics/blob/master/tutorials/treePL.md)
	- General pipeline for performing phylogenetic divergence dating with TreePL



## Scripts
These are minor scripts that may be useful, but don't warrant their own GitHub page or recognition as a full pipeline. These can be found within this repository. Most of the dependencies can be (and I would recommend that they are) installed with [Anaconda](https://www.anaconda.com/products/individual)/[Mamba](https://mamba.readthedocs.io/en/latest/) or with CRAN in [R](https://cran.r-project.org/).

- **`bgb_to_treedata.R`**: Modified script to convert [BioGeoBEARS](http://phylo.wikidot.com/biogeobears) results into [treeio](https://yulab-smu.top/treedata-book/)- and [ggtree](https://yulab-smu.top/treedata-book/)-readable treedata. Original script by [RevGadgets/Cody Howard/Carrie M Tribble](https://github.com/revbayes/RevGadgets/issues/9) 
	- Dependencies: [BioGeoBEARS](http://phylo.wikidot.com/biogeobears), [treeio](https://yulab-smu.top/treedata-book/), [RevGadgets](https://revbayes.github.io/tutorials/intro/revgadgets), [naniar](https://naniar.njtierney.com/)
- **`ConcatAln.py`**: Modified Python script to concatenate a directory of gene alignments for phylogenetics. Original script from [Andrew Mason](https://github.com/masonaj157)
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages), [tqdm](https://github.com/tqdm/tqdm)
- **`ConcatSNPs.py`**: Script will read a VCF and convert it to a concatenated fasta and nexus for input into SVDQuartets. It will produce separate sequences for each allele (i.e., Sample1_a1, Sample1_a2) and will include a partition to coalesce the two alleles together in SVDQuartets
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages), [PyVCF](https://pyvcf.readthedocs.io/en/latest/INTRO.html)
- **`create_hyde_triples.R`**: R function to create a file of all possible triples for [HyDe](https://github.com/pblischak/HyDe). *triple* = parent_1 --> hybrid <-- parent_2
	- Dependencies: NA
- **`create_pairwise_combos.R`**: R function to create a file of all possible pairwise comparisons given vector of options.
	- Dependencies: NA
- **`FastaRenamer.py`**: Python script to rename fasta files sequentially. Blatantly plagarized from [Andrew Mason](https://github.com/masonaj157).
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`fqdump.py`**: Python script that acts as a wrapper for fastq-dump or fasterq-dump in order to get around their dumb default settings and to convert all quality scores to PHRED+33
	- [Manual](https://github.com/RhettRautsaw/Bioinformatics/blob/master/scripts/script_docs/fqdump.md)
	- Depdendencies: [Python](https://www.python.org/), [sra-tools](https://github.com/ncbi/sra-tools), [pigz](https://zlib.net/pigz/), [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
- **`gif_maker.py`**: Python script which converts mov files to gifs. This is useful for converting screen-recordings (done with macOS Quicktime) into gifs for GitHub pages (see my [QGIS Tutorial](https://github.com/RhettRautsaw/GIS_Tutorial/blob/master/QGIS_Tutorial.md))
	- Dependencies: [Python](https://www.python.org/), [ffmpeg](https://ffmpeg.org/), [gifsicle](https://github.com/kohler/gifsicle)
- **`Linearize.py`**: Python script which takes multiline fasta and makes it a single-line fasta and vice-versa.
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`MCL.py`**: Script to run a self-BLAST search and cluster sequences using a Markov Clustering. Currently designed to keep only the sequences in the largest cluster. **Need to modify to allow users to choose to keep all sequences or only largest cluster**
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [MCL](https://micans.org/mcl/)
- **`PresAbsCheck.py`**: Python script to check if transcripts are "present" based on read coverage across the transcript. Sequences with low coverage across a certain percentage of the transcript may be indicative of a misassemply or false-positive annotation. Default is to check for transcripts with <5x coverage for >10% of the total transcript length. 
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages), [Numpy](https://numpy.org/), [Pandas](https://pandas.pydata.org/), [dfply](https://github.com/kieferk/dfply), [bwa](https://github.com/lh3/bwa), [samtools](http://www.htslib.org/), [bedtools](https://bedtools.readthedocs.io/en/latest/).
- **`RemDupRemAmb.py`**: Python script to remove sequences with ambiguities or duplicates. Blatantly plagarized from [Andrew Mason](https://github.com/masonaj157).
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`RemoveDups.py`**: Python script which takes a fasta file and removes sequences with duplicate names (keeping the first sequence it encounters).
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`RemoveSeqs.py`**: Python script which takes a fasta file and list of sequences and removes these sequences from the fasta.
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`RemoveStop.py`**: Python script which takes a transcriptome fasta and removes the last three base pairs (*i.e.*, stop codon).
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`RenameDups.py`**: Python script which takes a fasta file and searches for duplicate sequence IDs, adding a unique index if found.
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`Rogue_CDS_Finder.py`**: Python script designed help detect highly expressed ORFs/CDSs that may have been missed during annotation or accidentally removed during another processing step (*i.e.*, Rogue CDSs). Blatantly plagarized from [Andrew Mason](https://github.com/masonaj157).
	- [Manual](https://github.com/RhettRautsaw/Bioinformatics/blob/master/scripts/Rogue_CDS_Finder)
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- **`RSEM_combiner.R`**: Unix-executable R script to combine RSEM results for multiple individuals into a single dataframe/csv
	- Dependencies: [R argparse](https://cran.r-project.org/web/packages/argparse/vignettes/argparse.html), [readr](https://readr.tidyverse.org/), [stringr](https://stringr.tidyverse.org/), [tidyr](https://tidyr.tidyverse.org/), [dplyr](https://dplyr.tidyverse.org/), [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
- **`SCP_Deliver.py`**: This is the first Python script I ever wrote. All it does is remember that pesky `scp` formatting for you. Not really that useful, but it's nostalgic for me.
	- Dependencies: NA
- **`SortFasta.py`**: Python script to sort a fasta by either sequence names or length. 
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`Ungap.py`**: Python script to remove gaps from a fasta file/alignment.
	- Dependencies: [Biopython](https://biopython.org/wiki/Packages)
- **`VCF_subsetR.R`**: Unix-executable R script designed to take a vcf/bcf file and subset it according to provided metadata

<br><br>
***
<br><br>

# Additional Resources
## Misc Functions
These are functions that I have found useful over the years, but seem like they would be really easy to forget. So I've included them here. You might find them useful as well, but they are mostly here for me. 

**UNIX**
- Parallel Compress Directory:
	- `tar -cvf - directory | pigz -9 -p 20 > directory_archive.tar.gz`
- Recursively Find File by Name and Remove:
	- `find . -name "name*" -exec rm {} +`
- List file hosted on remote https
	- `lftp -u username,password\!\! https://www.website.com`
	- `du -a > manifest.txt`
- Entrez Direct e-utils grab accession information
	- `esearch -db sra -query SRR12915634 | esummary | xtract -pattern DocumentSummary -element Run@acc Experiment@acc Sample@acc Biosample Bioproject`
- Allow copy-pasting multiple lines in Python
	- `echo "set enable-bracketed-paste off" >> ~/.inputrc`

**REGEX**
- Remove duplicated text (filter for unique)
	- `perl -pe 's/\b(\w+)(?:,\s+\1\b)+/$1/g'`: Replace Duplicated Words when words are back to back. 
	- `perl -pe 's/\b(\w+)\b\s*,\s*(?=.*\1)//g'`: Replace Duplicated words per line`

**Excel/Google Sheets**
- VLOOKUP
	- `=VLOOKUP(CELL, RANGE, COLUMN, FALSE)`
- VLOOKUP concatenate multiple matches
	- `ArrayFormula(TEXTJOIN("; ", TRUE, IF(A1=AnotherSheet!B:B, D:D, "")))`

***

## Links
These are websites/links that I've found useful at various times through the years. 
- [Rhett's Starred GitHub Repositories](https://github.com/RhettRautsaw?tab=stars)
- [Filtering SNPs](https://otagomohio.github.io/2019-06-11_GBS_EE/sessions/filteringSNPs.html)
- [SNP Filtering Tutorial](http://www.ddocent.com/filtering/)
- [Phyluce: Loci Harvesting](https://phyluce.readthedocs.io/en/latest/tutorial-four.html)
- [PLWS: Phylogenomics from Low-coverage Whole-genome Sequencing](https://github.com/xtmtd/PLWS)
- [Guide to Tmux](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)
- https://www.robertlanfear.com/blog/files/short_read_mappers.html

