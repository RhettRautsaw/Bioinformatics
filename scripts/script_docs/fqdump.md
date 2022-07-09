# fqdump
## Rhett M. Rautsaw

***

`fastq-dump` and `fasterq-dump` should really be called **fastq-dum(b)** because they have dumb defaults that screw with downstream processing. Use my script `fqdump.py` and let me help format those sequences and change the quality scores for you into something that makes sense and won't screw with downstream data processing. 

## Installation
```
conda install parallel sra-tools bbmap pigz
```

## Use
```
fqdump.py -s SRX6081244 -n Crhod-SRS4983010_RNA_VG -t 16
```
So EASY!!!

<br>

*** 
<br>

## fastq-dump

If you don't want to use my script or for some reason it is not working...fine...but at least read what I have here about the settings for `fastq-dump`.

This is my recommended settings:
```
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --defline-seq '@$ac:$si $ri' --defline-qual '+' SRX
```

### Rationale
#### **Defline Seq**
The `--defline-seq` is the important part. So lets see what the fastq looks like when you leave that off.
```
fastq-dump --gzip --skip-technical --readids --read-filter pass  --dumpbase --split-3 --clip SRX6081244 | head -4

@SRR9313681.1.1 FCH2N7KBBXX:7:1101:5171:1457 length=100
TTGTTGCATTAAATGTGCTTTGCNATTGGCAAATGGTAAAATTGCCATTTTTTTTTCTTGACAACAAATAAGCCATTGTTGTTGTACGTGCCCATGGAGC
+
``eeejjjjejjejjjjjejeejBjjjjjjVejjjjjjjjjjjjjjjjjjjjjjjjKK`KVejj`ejjjjejKKKKV[ej`jjej``[ejVG``jVe`e`
```

I don't know why this is the default, but you don't need most of that information and the periods cause problems for several downstream applications like `Trim-Galore` and `Trinity`. You can use `--defline-seq` to fix this and get your custom fastq header output. So what are your options for a fastq header?

| variable | description          |
|----------|----------------------|
| $ac      | accession            |
| $si      | spot id              |
| $sn      | spot name            |
| $sg      | spot group (barcode) |
| $sl      | spot length in bases |
| $ri      | read number          |
| $rn      | read name            |
| $rl      | read length in bases |

Let's check them all out

```
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --defline-seq '$ac | SI = $si | SN = $sn | SG = $sg | SL = $sl | RI = $ri | RN = $rn | RL = $rl' --defline-qual '+' SRX6081244 | head -4

@SRR9313681 | SI = 1 | SN = FCH2N7KBBXX:7:1101:5171:1457 | SG = ATCACGAT | SL = 200 | RI = 1 | RN = forward | RL = 100
TTGTTGCATTAAATGTGCTTTGCNATTGGCAAATGGTAAAATTGCCATTTTTTTTTCTTGACAACAAATAAGCCATTGTTGTTGTACGTGCCCATGGAGC
+
``eeejjjjejjejjjjjejeejBjjjjjjVejjjjjjjjjjjjjjjjjjjjjjjjKK`KVejj`ejjjjejKKKKV[ej`jjej``[ejVG``jVe`e`
```

So really all we need is the accession number (`$ac`), a unique identifier (`$si`), and the readid (`$ri`). And for the quality header, we can just keep it simpler with '+'.
```
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --defline-seq '@$ac:$si $ri' --defline-qual '+' SRX6081244 | head -4

@SRR9313681:1 1
TTGTTGCATTAAATGTGCTTTGCNATTGGCAAATGGTAAAATTGCCATTTTTTTTTCTTGACAACAAATAAGCCATTGTTGTTGTACGTGCCCATGGAGC
+
``eeejjjjejjejjjjjejeejBjjjjjjVejjjjjjjjjjjjjjjjjjjjjjjjKK`KVejj`ejjjjejKKKKV[ej`jjej``[ejVG``jVe`e`
```
#### **Quality Scores**
The other problem is that sometimes older datasets can be in `phred64` format instead of `phred33` format. This can similarly cause problems for trimming sequences. You can tell when a dataset is in `phred64` by the presence of lowercase letters in the quality scores.

We can use `bbmap`'s `reformat.sh` tool to convert quality scored into `phred33` format. 
```
reformat.sh in=SRR9313681_pass_1.fastq.gz out=SRR9313681_pass_1_reformat.fastq.gz qout=33 overwrite=true
```