# Bioinformatics
## Rhett M. Rautsaw
This repository contains many simple, but potentially useful bioinformatics scripts as well as functions and resources for improving your skillset and performing really cool analyses.

## Scripts

- **`create_hyde_triples.R`**: R function to create a file of all possible triples for [HyDe](https://github.com/pblischak/HyDe). *triple* = parent_1 --> hybrid <-- parent_2
- **`create_pairwise_combos.R`**: R function to create a file of all possible pairwise comparisons given vector of options.
- **`Linearize.py`**: Python script which takes multiline fasta and makes it a single-line fasta.
- **`RemoveDups.py`**: Python script which takes a fasta file and removes sequences with duplicate names (keeping the first sequence it encounters).
- **`RemoveStop.py`**: Python script which takes a transcriptome fasta and removes the last three base pairs (*i.e.*, stop codon).
- **`RenameDups.py`**: Python script which takes a fasta file and searches for duplicate sequence IDs, adding a unique index if found.

## Functions

### UNIX
- Parallel Compress Directory:
	- `tar -cvf - directory | pigz -9 -p 20 > directory_archive.tar.gz`
- Recursively Find File by Name and Remove:
	- `find . -name "name*" -exec rm {} +`
- List file hosted on remote https
	- `lftp -u username,password\!\! https://www.website.com`
	- `du -a > manifest.txt`

### Excel/Google Sheets
- VLOOKUP concatenate multiple matches
	- ArrayFormula(TEXTJOIN("; ", TRUE, IF(A1=AnotherSheet!B:B, D:D, "")))

## Scripting Resources

- [Scripting OS X](https://scriptingosx.com/witchcraft/)
	- [Moving to zsh](https://scriptingosx.com/2019/06/moving-to-zsh/)
	- [OhMyZsh](https://ohmyz.sh/)

## Analysis Resources

- [Filtering SNPs](https://otagomohio.github.io/2019-06-11_GBS_EE/sessions/filteringSNPs.html)
- [SNP Filtering Tutorial](http://www.ddocent.com/filtering/)
- [Phyluce: Loci Harvesting](https://phyluce.readthedocs.io/en/latest/tutorial-four.html)
- [PLWS: Phylogenomics from Low-coverage Whole-genome Sequencing](https://github.com/xtmtd/PLWS)
- [Guide to Tmux](https://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)
- [ENMeval](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12261)
	- [Vignette](https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html)
- [kuenm](https://peerj.com/articles/6281/)
