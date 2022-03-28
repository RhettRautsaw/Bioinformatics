# HiFi Genome Assembly and Annotation
## Rhett M. Rautsaw

***

This markdown document walks through the general pipeline followed for [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) genome assembly and annotation. I use our genome for *Crotalus lepidus* (Rock Rattlesnake; Clepi-CLP2201_WGS_blood) as an example. 

# Quality Check and Stats
We can use [NanoPlot](https://github.com/wdecoster/NanoPlot) to check out some general statisitics on our long-read data such as the average/median read lengths and N50 before assembly.
```
NanoPlot -t 80 -o NanoPlot_res --fastq Clepi-CLP2201_WGS_blood_hifi1.fastq.gz Clepi-CLP2201_WGS_blood_hifi2.fastq.gz
```

# Assembly
After testing several methodologies including (*i.e.*, [HiCanu](https://github.com/marbl/canu), [Peregrine](https://github.com/cschin/peregrine-2021), [MaSuRCA](https://github.com/alekseyzimin/masurca), and [hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html)), we chose to use [hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html) as our primary assembler as it was both the fastest and produced the best assembly. 
```
hifiasm -o Clepi-CLP2201_WGS_blood -t 80 Clepi-CLP2201_WGS_blood_hifi1.fastq.gz Clepi-CLP2201_WGS_blood_hifi2.fastq.gz
```

hifiasm outputs the primary assembly in `gfa` format, we can convert this to fasta format using `awk`

```
awk '/^S/{print ">"$2;print $3}' Clepi-CLP2201_WGS_blood.bp.p_ctg.gfa > Clepi-CLP2201_WGS_blood.p_ctg.fasta
```

hifiasm outputs many different files, you can check out what all the different files are [here](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html).
Now we can check out how our assembly did using [bbstats](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
```
bbstats -in=Clepi-CLP2201_WGS_blood.p_ctg.fasta &> Clepi-CLP2201_WGS_blood.p_ctg.fasta.stats.txt
```

# BUSCO
To assess the quality of your genome assembly, you will want to run [BUSCO](https://busco.ezlab.org/). Ideally, you will have BUSCO completeness at least >90%
```
busco -i Clepi-CLP2201_WGS_blood.p_ctg.fasta -m genome -l /PATH/TO/BUSCO/DB/2020-09-10_tetrapoda_odb10 -c 80 -o 06_busco
```

# Annotation
## TE & Repeat Annotation/Masking
To identify and annotate transposable elements, long-terminal repeats (LTRs), and other repeats in the genome, we use the [Extensive De-novo TE Annotator (EDTA)](https://github.com/oushujun/EDTA) which performs [RepeatModeler/RepeatMasker]() within it if `--sensitive 1` is turned on. **Warning: `--sensitive 1` is slow and will take a while!**
```
perl /zfs/venom/Rhett/bin/EDTA/EDTA.pl --genome ../04_hifiasm/Clepi-CLP2201_WGS_blood.p_ctg.fasta --threads 80 --sensitive 1 --anno 1
```

## Soft-masking
[EDTA](https://github.com/oushujun/EDTA) automatically hard-masks your genome – meaning that all repeat elements get converted from nucleotides to "N". However, most annotation programs prefer that your genome is soft-masked – meaning that capital letters get converted to lowercase letters where repeat elements exist. We can use the output of [EDTA](https://github.com/oushujun/EDTA) to perform a soft-mask.

This is a multi-step process. First rename your hard-masked fasta so that you can identify as the hard-masked version.
```
mv Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.masked Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.hard.masked
```

Now create a new directory where we can do the masking and copy in the genome and relevant masking elements
```
mkdir softmask 
cp ../04_hifiasm/Clepi-CLP2201_WGS_blood.p_ctg.fasta softmask/
cp Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.anno/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out softmask/
cd softmask
```

Now we are ready to actually do the soft-masking; mimicking how [EDTA](https://github.com/oushujun/EDTA) did the hard mask.
```
perl ~/.conda/envs/EDTA/share/EDTA/util/make_masked.pl -genome Clepi-CLP2201_WGS_blood.p_ctg.fasta -rmout Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0 -misschar N -threads 80

mv Clepi-CLP2201_WGS_blood.p_ctg.fasta.new.masked ../Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked

rm Clepi-CLP2201_WGS_blood.p_ctg.fasta Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out
```

## LAI score
We can also use the [EDTA](https://github.com/oushujun/EDTA) output to calculate the [LTR-Assembly Index (LAI)](https://academic.oup.com/nar/article/46/21/e126/5068908) which is another genome quality index. 
<center>

| Category  | LAI   |
|-----------|-------|
| Draft     | < 10  |
| Reference | 10-20 |
| Gold      | ≥ 20  |

</center>

```
LAI -genome Clepi-CLP2201_WGS_blood.p_ctg.fasta -intact Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.raw/LTR/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.pass.list -all Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.anno/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.out -t 80
```


## Gene Annotation
To do gene annotation, we can use [funannotate](https://funannotate.readthedocs.io/en/latest/#) which is automated gene prediction and functional annotation package. This program is very well documented, so you can modify for your purposes; however, we will generally follow 4 easy-to-use steps. 
1. Training (using RNA data)
2. Prediction (using training data and additional protein databases)
3. Updating (using RNA data to back over predictions and add UTR annotations)
4. Functional Annotation (annotate the final predictions with actual gene names, GO Terms, and functional information)

### Installation
I've found that installation can be a little tricky, but far easier than most other genome annotation programs. I recommend using Mamba to speed things up
```
conda create -n funannotate mamba
conda  activate  funannotate
mamba install "python>=3.6,<3.9" funannotate
mamba install -c bioconda perl-bioperl
mamba install -c bioconda eggnog-mapper

# Create folder in which we can install databases for funannotate
mkdir /MODIFY/THIS/PATH/bin/funannotate
cd /MODIFY/THIS/PATH/bin/funannotate
funannotate setup -d funannotate_db -b tetrapoda

# There appears to be a folder missing in the conda environment which is required for EggNogg. Create this directory and then download the necessary EggNog database
mkdir /home/rrautsa/.conda/envs/funannotate/lib/python3.8/site-packages/data
download_eggnog_data.py

# Lastly, InterProScan5 is not available via conda. I'd recommend installing this locally. This will take a while though to download. 
mkdir interproscan
cd interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.54-87.0/interproscan-5.54-87.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.54-87.0/interproscan-5.54-87.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.54-87.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.54-87.0-*-bit.tar.gz
python3 initial_setup.py

# You can check that everything is install properly with the following commands. The test will take a while.
funannotate check --show-versions
funannotate test -t all --cpus 16
```

### 1. Training
```
funannotate train -i ../07_repeat/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked -o annotate \
	--left ../00_rna/Clepi-CLP1932_RNA_heart_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_liver_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_muscle_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_VG_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_kidney_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_pancreas_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_VG_R1_trim.fastq.gz \
	--right ../00_rna/Clepi-CLP1932_RNA_heart_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_liver_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_muscle_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_VG_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_kidney_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_pancreas_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_VG_R2_trim.fastq.gz \
	--no_trimmomatic --max_intronlen 30000 \
	--cpus 80 --species "Crotalus lepidus"
```

### 2. Prediction

```
funannotate predict -i ../07_repeat/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked -o annotate \
	--transcript_evidence ../00_rna/Clepi_concensus_transcriptome_97_v1_renamed.fasta \
	--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
	--busco_db tetrapoda --busco_seed_species chicken --organism other --max_intronlen 30000 \
	--cpus 80 --species "Crotalus lepidus"
```

### 3. Updating

```
funannotate update -i annotate --cpus 80
```

### 4. Functional Annotation

```
funannotate iprscan -i annotate -c 80 -m local --iprscan_path /zfs/venom/Rhett/bin/funannotate/interproscan/interproscan.sh
funannotate annotate -i annotate --busco_db tetrapoda --cpus 80
```


# Full PBS Script
```
#PBS -N hifi_processing
#PBS -l select=1:ncpus=80:mem=1500gb,walltime=336:00:00 
#PBS -q venom
#PBS -j oe
#PBS -m abe
#PBS -M rrautsa@g.clemson.edu

module load anaconda3/2019.10-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

# QUALITY CHECK & STATS
NanoPlot -t 80 -o NanoPlot_res --fastq 2021-10-22_UDel_PacBio_HiFi-1/01_concat/Clepi-CLP2201_WGS_blood_hifi.fastq.gz 2021-11-12_UDel_PacBio_HiFi-2/01_concat/Clepi-CLP2201_WGS_blood_hifi.fastq.gz

# 04 ASSEMBLE
mkdir 04_hifiasm
cd 04_hifiasm
hifiasm -o Clepi-CLP2201_WGS_blood -t 80 2021-10-22_UDel_PacBio_HiFi-1/01_concat/Clepi-CLP2201_WGS_blood_hifi.fastq.gz 2021-11-12_UDel_PacBio_HiFi-2/01_concat/Clepi-CLP2201_WGS_blood_hifi.fastq.gz
awk '/^S/{print ">"$2;print $3}' Clepi-CLP2201_WGS_blood.bp.p_ctg.gfa > Clepi-CLP2201_WGS_blood.p_ctg.fasta
bbstats -in=Clepi-CLP2201_WGS_blood.p_ctg.fasta &> Clepi-CLP2201_WGS_blood.p_ctg.fa.stats.txt
source deactivate
cd ..

# 06 BUSCO
source activate busco521_env
busco -i 04_hifiasm/Clepi-CLP2201_WGS_blood.p_ctg.fasta -m genome -l /zfs/venom/02_databases/busco/2020-09-10_tetrapoda_odb10 -c 80 -o 06_busco
source deactivate

# 07 ANNOTATION
## EDTA: TE ANNOTATION & REPEAT MASKING
source activate EDTA
mkdir 07_repeat
cd 07_repeat
perl /zfs/venom/Rhett/bin/EDTA/EDTA.pl --genome ../04_hifiasm/Clepi-CLP2201_WGS_blood.p_ctg.fasta --threads 80 --sensitive 1 --anno 1

## SOFT MASK
mv Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.masked Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.hard.masked
mkdir softmask 
cp ../04_hifiasm/Clepi-CLP2201_WGS_blood.p_ctg.fasta softmask/
cp Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.anno/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out softmask/
cd softmask
perl ~/.conda/envs/EDTA/share/EDTA/util/make_masked.pl -genome Clepi-CLP2201_WGS_blood.p_ctg.fasta -rmout Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0 -misschar N -threads 80
mv Clepi-CLP2201_WGS_blood.p_ctg.fasta.new.masked ../Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked
rm Clepi-CLP2201_WGS_blood.p_ctg.fasta Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.RM.out
cd ..

## LAI
mkdir LAI; cd LAI
LAI -genome ../Clepi-CLP2201_WGS_blood.p_ctg.fasta -intact ../Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.raw/LTR/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.pass.list -all ../Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.EDTA.anno/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.out -t 80
cd ../..
source deactivate

## FUNANNOTATE
module load anaconda3/2019.10-gcc/8.3.1
source activate funannotate
mkdir 07_funannotate; cd 07_funannotate

### TRAIN
funannotate train -i ../07_repeat/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked -o annotate \
	--left ../00_rna/Clepi-CLP1932_RNA_heart_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_liver_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_muscle_R1_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_VG_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_kidney_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_pancreas_R1_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_VG_R1_trim.fastq.gz \
	--right ../00_rna/Clepi-CLP1932_RNA_heart_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_liver_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_muscle_R2_trim.fastq.gz ../00_rna/Clepi-CLP1932_RNA_VG_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_kidney_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_pancreas_R2_trim.fastq.gz ../00_rna/Clepi-CLP2201_RNA_VG_R2_trim.fastq.gz \
	--no_trimmomatic --max_intronlen 30000 \
	--cpus 80 --species "Crotalus lepidus"

### PREDICT GENES
funannotate predict -i ../07_repeat/Clepi-CLP2201_WGS_blood.p_ctg.fasta.mod.MAKER.soft.masked -o annotate \
	--transcript_evidence ../00_rna/Clepi_concensus_transcriptome_97_v1_renamed.fasta \
	--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
	--busco_db tetrapoda --busco_seed_species chicken --organism other --max_intronlen 30000 \
	--cpus 80 --species "Crotalus lepidus"

### UPDATE UTRs
funannotate update -i annotate --cpus 80

### FUNCTIONAL ANNOTATION
funannotate iprscan -i annotate -c 80 -m local --iprscan_path /zfs/venom/Rhett/bin/funannotate/interproscan/interproscan.sh
funannotate annotate -i annotate --cpus 80
```
