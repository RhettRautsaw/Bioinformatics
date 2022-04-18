# How to use Palmetto Cliffnotes

## Login
Logging into Palmetto is as easy as:
```
ssh USERNAME@login.palmetto.clemson.edu
```

> Note that when you login in you are legitimately on what is called the `login` node. Absolutely no computation other than navigating folders, writing scripts, and submitting jobs should be done from the `login` node. To do actual computation, you should start an interactive job or write a script and submit it. 

## Important Directories
- `/home/USERNAME` or `~` 
	- Home Directory. This is where your conda environments live as well as your bashrc
	- You can also install things in this directory and can hold data here, but you have limited storage space (100gb), so use it wisely. 
	- This directory is also backed up regularly.
- `/scratch1/USERNAME`
	- Scratch Directory. This is where analyses and everything should be done.
	- Generally it is best to transfer your data here with **Globus** or `scp`, run your analyses, and then transfer back. If you need to transfer files, check out this page on how to use `scp` and/or **Globus**. https://www.palmetto.clemson.edu/palmetto/basic/started/
	- This directory is not backed up, so if you lose something in here, it is gone forever.

## PBS Scripts
To write scripts (which you can do in the `login` node), I use the command `nano` which is a command-line based text editor. Simply use `nano JOBNAME.pbs` to create and start writing the file. Use `Ctrl + O` to save the file and `Ctrl + X` to exit the file. 

As shown above, jobs are kept in `pbs` scripts which start with the headers shown below:
```
#PBS -N IQTREE
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=fdr,walltime=72:00:00
#PBS -j oe
#PBS -M USERNAME@g.clemson.edu
#PBS -m abe
```
What does all this mean? Well it's a set of instructions for what resources Palmetto's `scheduler` should grab for you and specifications for how to output the results and inform you about them. 
- -N = Job name
- -l = Job specifications
	- `select`: The number of nodes/computers you want to use. You don't need more than one unless you know how to parallelize across nodes (which I will show you later)
	- `ncpus`: The number of cores on each select node
	- `mem`: The amount of memory on each select node
	- `interconnect`: Palmetto is separated into what they call "clusters". There are three main options here:
		- `1g`: This is the phase with the smallest nodes (~ 8 cpus & 15-30 gb memory). However there is one small group of computers in `1g` that is 16 cpus and 251 gb memory. This group in particular can be really useful. 
		- `fdr`: This is the phase with the medium-sized nodes (~ 16-28 cpus & 62-125gb memory).
		- `hdr`: This is the phase with medium-large sized nodes (~ 40 cpus & 372gb memory)
	- `walltime`: This is the amount of time your job will take. Although the Palmetto people don't like this, I generally just always put this at the maximum time which is 72:00:00. 
- -j = This one is optional, but this joins the resulting output and error of the job into one file. 
- -M = Email address to send messages to.
- -m = What type of messages should the job send you
	- a = abort
	- b = begin
	- e = end

If your job needs more memory or time, you need to use the `bigmem` nodes (`walltime: 168:00:00`). You can specify the `bigmem` nodes with the `-q` command. Note that if you do this, then you do not need to specify `interconnect` in the `-l` command.
```
#PBS -N IQTREE
#PBS -l select=1:ncpus=16:mem=62gb,walltime=72:00:00
#PBS -q bigmem
#PBS -j oe
#PBS -M USERNAME@g.clemson.edu
#PBS -m abe
```

## Optimizing Resources
When submitting a job, assuming you are not extremely picky about your computational resources, it is best to use the command `whatsfree` to check what resources are available/FREE and then modify you PBS job accordingly. This will ensure that your job starts soon after you submit it. 

## Submitting/Monitoring Jobs

### Submit Job
```
qsub JOB_FILE.pbs
```

### Interactive Jobs
It is recommended that you write `pbs` scripts to keep a log of what you have run. But you can also run things interactively by using `-I` flag. Below you will see that instead of specifying a `pbs` script file, I specified the `-l` command which matches the PBS script above.
```
qsub -I -l select=1:ncpus=16:mem=62gb:interconnect=fdr,walltime=72:00:00
``` 

### Monitor Jobs
After you submit a job, you can monitor it with the following commands.

Check jobs submitted under your username
```
qstat -u USERNAME
```

Peek at what is happening inside that job. 
```
qpeek <job id number>
```

### Kill Job
If a job is not running properly and you know that, you can kill the job.
```
qdel <job id number>
```

## Modules
Palmetto has many pre-installed modules that you can load. You can check them out using
```
module avail
```

And load them into your PATH for use with
```
module load anaconda3/2021.05-gcc/8.3.1
```
Above I loaded the anaconda module because anaconda makes package installation easy. 


## Put it All Together
To walk through an example, lets infer a phylogeny based on a concatenated alignment with gene partitions. First transfer your alignment and partition file to the `scratch1/USERNAME` directory. Then login to Palmetto and navigate to that directory.

When you login, you are on the `login` node. Start an interactive job to create an anaconda environment with iqtree in it.
```
qsub -I -l select=1:ncpus=16:mem=62gb:interconnect=fdr,walltime=72:00:00
```

You should be "teleported" to a different computer and should see a change in your terminal prompt from `login` to `node####`. But what you also may notice is that when you "teleported", you also got moved back to your HOME directory. This is important to realize, because **every time** you submit a job, you are being "teleported" to a new node and back to your HOME directory. Luckily, Palmetto remembers where you started from and it is easy to navigate back there by calling the variable `$PBS_O_WORKDIR`. 
 
Now you can load anaconda and create the iqtree environment
```
module load anaconda3/2021.05-gcc/8.3.1
conda create -n iqtree_env iqtree biopython parallel
```

We can also write our script while in an interactive job. Use nano to write the following script:
```
#PBS -N IQTREE
#PBS -l select=1:ncpus=24:mem=1000gb,walltime=72:00:00
#PBS -q bigmem
#PBS -j oe
#PBS -M USERNAME@g.clemson.edu
#PBS -m abe

module load anaconda3/2021.05-gcc/8.3.1
conda activate iqtree_env 
# you may need to use this instead:
# source activate iqtree_env

cd $PBS_O_WORKDIR

iqtree -s SEQUENCE_FILE.fasta -T 24 -B 1000 -p PARTITION_FILE.txt
```

Notice in the script that I first loaded anaconda and activated my conda environment. Because the script is logging you into a different node/computer, you also lose anything you may have loaded in your interactive job or on the `login` node. I then similarly navigated back to the directory in which my script is held with `cd $PBS_O_WORKDIR`. 

Finally, I can then run iqtree on that computational node. 

## GNU-Parallel
Lets say you have 1000 genes that you want to infer phylogenies for separately. Lets also say that each of these genes has their own partitions on which you want to allow different evolutionary models to be fit. Seems like it could be difficult to code right? Wrong...

This can easily be done with the help of GNU-Parallel. First create a tab-delimited list of all your genes/paritions. Genes in column 1 and Partitions in column 2:
```
ls *.fasta > GENES.txt
ls *.txt > PARTITIONS.txt
paste GENES.txt PARTITIONS.txt > FINAL_LIST.txt
```

You can then either specify to parallelize across multiple cpus of a single node:
```
#PBS -N IQTREE
#PBS -l select=1:ncpus=24:mem=1000gb,walltime=72:00:00
#PBS -q bigmem
#PBS -j oe
#PBS -M USERNAME@g.clemson.edu
#PBS -m abe

module load anaconda3/2021.05-gcc/8.3.1
conda activate iqtree_env 

cd $PBS_O_WORKDIR

parallel -a FINAL_LIST.txt -j 24 --colsep '\t' --bar '
	iqtree -s {1} -T 1 -B 1000 -p {2}
	'
```

Or across multiple nodes (recommended):
```
#PBS -N IQTREE
#PBS -l select=200:ncpus=8:mem=30gb:interconnect=1g,walltime=72:00:00
#PBS -j oe
#PBS -M USERNAME@g.clemson.edu
#PBS -m abe

module load anaconda3/2021.05-gcc/8.3.1
conda activate iqtree_env 

cd $PBS_O_WORKDIR

parallel -a FINAL_LIST.txt -j 1 --colsep '\t' --bar --workdir $PWD --sshloginfile $PBS_NODEFILE '
	module load anaconda3/2021.05-gcc/8.3.1
	conda activate iqtree_env 
	iqtree -s {1} -T 8 -B 1000 -p {2}
	'
```
Some commonalities in these two examples is that I called `parallel` and specified my tab-delimited list with the `-a` flag and `--colsep '\t'`. I also added a progress bar with the `--bar` flag. I also specified the insertion of the column 1 values vs column 2 values in my iqtree command with `{1}` and `{2}` in both. 

Now some differences...

In the first example, I specified `-j 24` in parallel and reduced the number of threads `-T` in IQTree. Now each gene will only get 1 thread, but I can run 24 genes at one time. Depending on the process being done and efficiency of the program, this may even out in terms of time at the end of the day, but it is often faster. However, you have to be cautious regarding memory usage when using this technique. 

Because of this, I generally recommend the latter technique of parallelizing across multiple nodes/computers each with it's own dedicated memory. In this second example, I asked for many smaller nodes instead of one larger node. The first thing to note is that in addition to saving the variable `$PBS_O_WORKDIR`, Palmetto also saves `$PBS_NODEFILE` which is a list of all the computers available to you from the `-l select:##` flag in your PBS script. Parallel can use this list to `ssh` into each of the nodes and run jobs. In this example, because parallel is logging into 200 different nodes at one time, we need to reload our module and either navigate back to our working directory `cd $PBS_O_WORKDIR` or specify the `--workdir` flag to tell parallel to automatically navigate back to this folder. However, with this technique, I can leave my IQTree threads higher `-T 8` since each node only has a single job `-j 1`. Once one gene finishes, the next gene is automatically launched on the available node.

### Combination Parallelization
The last thing I will show you is that not only can you match up multiple parameters in a tab-separated list to be called together (as we did above). You can also specify multiple lists to do combinations of parameters. For example, lets say we wanted to run each of the 1000 genes with a different starting seed (random number starting point for ML optimization). I'm not sure why you'd want to do this, but it was the easiest thing I could think of as an example. Regardless, GNU-parallel also makes this easy. Simply specify two lists:

```
parallel -a GENES.txt -a SEEDS.txt -j 1 --bar '
	iqtree -s {1} -T 8 -B 1000 -seed {2}
	'
```

If your `SEED.txt` list was simply a file containing the values 1 -- 100, then parallel would run:
```
gene: 0001	seed: 001
gene: 0001	seed: 002
gene: 0001	seed: 003
...
gene: 1000	seed: 098
gene: 1000	seed: 099
gene: 1000	seed: 100
```