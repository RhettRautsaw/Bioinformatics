# Bioinformatic Scripts
A repository of potentially useful bioinformatics scripts.

**Scripts**:

- **`create_hyde_triples.R`**: R function to create a file of all possible triples for [HyDe](https://github.com/pblischak/HyDe). *triple* = parent_1 --> hybrid <-- parent_2
- **`create_pairwise_combos.R`**: R function to create a file of all possible pairwise comparisons given vector of options
- **`Linearize.py`**: Python script which takes multiline fasta and makes it a single-line fasta
- **`RemoveStop.py`**: Python script which takes a transcriptome fasta and removes the last three base pairs (*i.e.*, stop codon)
- **`RenameDups.py`**: Python script which takes a fasta file and searches for duplicate sequence IDs, adding a unique index if found.
