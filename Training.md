# Bioinformatic Training Resources/Advice
## Rhett M. Rautsaw

If you don't have much experience in **unix**, [**python**](https://www.python.org/), [**R**](https://www.r-project.org/) [**(Studio)**](https://rstudio.com/), or bioinformatics in general, then this document has some general resources and information that will help you get started.

# Table of Contents
- [Training Resources](#training-resources)
- [Lists and Loops](#lists-and-loops)
	- [GNU-Parallel](#gnu-parallel)
- [Piping & Regular Expressions](#piping-and-regular-expressions)
- [Renaming Files](#renaming-files)
- [UNIX Profiles](#unix-profiles)
- [Software Management](#software-management)

# Training Resources

[Software Carpentry](https://software-carpentry.org/lessons/) and [Data Carpentry](https://datacarpentry.org/lessons/) offer several freely-available lessons/workshops to get
you started learning **unix**, [**python**](https://www.python.org/), [**R**](https://www.r-project.org/), and much more:

- [Unix Lesson](http://swcarpentry.github.io/shell-novice/)
- [Python Lesson 1](https://swcarpentry.github.io/python-novice-inflammation/)
- [Python Lesson 2](http://swcarpentry.github.io/python-novice-gapminder/)
- [R Lesson 1](http://swcarpentry.github.io/r-novice-inflammation/)
- [R Lesson 2](http://swcarpentry.github.io/r-novice-gapminder/)
- [Git Lesson](https://swcarpentry.github.io/git-novice/)
	- Git is used to keep track of versions of your files. It is really great for coding and is the basis of GitHub. It's actually pretty easy to use, so check out the workshop!

They even have workshops specific to genomics and ecology:

- [Genomic Workshops](https://datacarpentry.org/lessons/#genomics-workshop)
- [Ecology Workshops](https://datacarpentry.org/lessons/#ecology-workshop)

Some other useful resources might be:

- [Markdown Guide/Cheatsheet](https://www.markdownguide.org/getting-started/)
- [SWIRL: Learn R in R](https://swirlstats.com/)
- [Learn-Bioinformatics Resources](https://github.com/czbiohub/learn-bioinformatics)
- [The macOS School of Terminal Witchcraft and Wizardry](https://scriptingosx.com/witchcraft/)

# Lists and Loops

**Lists** are your *best friend* in coding. A list – very simply – is a file containing all the names or identifiers that you will want to **loop** through or process in [**parallel**](https://www.gnu.org/software/parallel/) (*e.g.*, one line for each of your samples). You can create a simple text file in **unix** with the command `nano list.txt` and save the file using keyboard shortcut `ctrl+o` and exit with `ctrl+x`. A list will look something like this:

``` bash
sample_1
sample_2
...
sample_n
```

You can also create lists by simply saving the standard output (STDOUT) of a function to a file. For example:

``` bash
ls *.fastq > list.txt
```

Once you create a list, you can provide it to a `for loop` or to the [**`parallel`**](https://www.gnu.org/software/parallel/) command to process many samples the same way – and even simultaneously – without actually having to redo the same command over and over again. For example:

``` bash
# For Loops will "loop" through your list one-by-one 
# processing each the same way
for i in $(cat list.txt)
do echo ${i}
done

# GNU-Parallel or just parallel, will process items simultaneously.
# In this example, parallel will process 3 items at a time.
parallel -a list.txt -j 3 "echo {}"
```

## GNU-Parallel



# Piping and Regular Expressions

**Piping** is the process of taking the STDOUT of one command and feeding it into the standard input (STDIN) of another. This is done with the vertical bar character `|` and can be useful for doing several commands in a row. Below, you will see that I use piping to create 3 different lists based on the fastq files I have in a directory.

The first list I create is a list of the unique individuals which we will use throughout the pipeline. In this list, each individual has it’s own line. The second two lists remove new line characters and replace them with spaces so that all individuals are on a single line.

``` bash
ls *.fastq.gz | cat | perl -pi -e 's/_.*.fastq.gz//g' | uniq > list.txt
sed "s/$/_R1.fastq/g" list.txt | tr '\n' ' ' | sed '/^\s*$/d' > list2.txt
sed "s/$/_R2.fastq/g" list.txt | tr '\n' ' ' | sed '/^\s*$/d' > list3.txt
```

Based on those lines of code, hopefully you understand that I am taking the output of the first command and feeding it into subsequent commands. However, you may be asking…what are `perl -pi -e`, `sed`, and `tr`?

Hopefully you have learned about `grep` and it’s ability to search for specific text in a document. What you may not know is that `grep` stands for Global Regular Expression Print. `grep`, `perl -pi -e`, `sed`, `tr`, and `awk` all use something called **Regular Expressions** or **regex** to find and/or replace text. I’ve used 3 different find/replace methods to generate my lists, but they are all essentially doing the same thing. regex takes a lot of practice, but they are incredibly useful in **unix** scripting. So take some time to learn them!

Useful regex resources:

- [RegexOne Training Lessons](https://regexone.com/)
- [Learn Regex](https://github.com/ziishaned/learn-regex)
- [Regex Cheatsheet](https://www.rexegg.com/regex-quickstart.html)
- [Regex Tester and Live Editor](https://regexr.com/)
- [Awk Tester and Live Editor](https://awk.js.org)
- [sed Tester and Live Editor](https://sed.js.org)

Each of the different regex methods (*i.e.*, `perl -pi -e`, `sed`, `awk`, etc.) have own features and unique formats; therefore, you may have to change between them. However, you will pick your favorite (whichever you learn first) and use it most frequently. My personal preference and the option I am most familiar with is `perl -pi -e`. [Here's an explanation of perl -pi -e.](https://stackoverflow.com/questions/6302025/perl-flags-pe-pi-p-w-d-i-t)

## Renaming Files

Sometimes files might have something appended to their name that you don’t want, so you need to rename them in specific way. If it is only one file, a simple `mv` command will do the trick. However, when there are a lot of files, renaming all of them can be a pain. Below, I provide three options for removing an unwanted underscore in all my file names.

**Option 1**: The `loop` and `mv` strategy

Here, we will loop through all files matching a pattern and `mv` them into a new file name. However, that new file name will be three arguments inside `${}` and separated by `/`. First the variable name again (`i`), then what you want to find (`sample_`), and finally what you want to replace that find with (`sample`).

``` bash
for i in sample_*.fastq.gz
	do mv $i ${i/sample_/sample}
	done
```

**Option 2**: The `find` and `rename` strategy

This option is very similar to the previous where we are finding all files with a certain pattern and then executing (`-exec`) the `rename` function on each. The rename function takes what we want to find (`sample_`), what we want to replace that find with (`sample`), and then the name of the file we want to perform this on. Since we are executing from `find`, the name of the file is represented by `{}`. This needs to be ended with `\;` to tell unix that your command is complete.

``` bash
find . -name "sample_*" -exec rename "sample_" "sample" {}\;
```

**Option 3**: The shell script strategy

This strategy takes a little more manual work, but is useful when you have less regular patterns to change, where each individual gets a slightly different name (*e.g.*, adding the species code which is found in a spreadsheet somewhere). In this situation, let Microsoft Excel help you.

To do this, create a list of all your individuals that need renamed. You can then copy-paste your list into Microsoft Excel to use things like find-replace and `VLOOKUP` to add a column including the species code.

You can also create a column to copy repeat elements like `mv` efficiently across all your individuals. Concatenate all your columns together and move it into a text document like the example below. Once you create your text bash script, you can just run `sh script.sh` to do all the work for you.

``` bash
#!/bin/bash

mv sample_001.fastq.gz sample-001.fastq.gz
mv sample-2.fastq.gz sample-002.fastq.gz
mv sample_3.fastq.gz sample-003.fastq.gz
#mv [MORE INDIVIDUALS]
```

## UNIX Profiles

Each step of the Guide assumes that you have already installed all the softwares that you need and it can be easily run by just typing in the appropriate command. Lets talk about how to make that assumption true.

First, it is important to know that unix has different **Shells**, also known as **command-line interpreters**. You've likely been using a Shell this whole time without knowing it because they are necessary to interpret every command you give. The standard for unix is the Bourne Shell (**sh**) or the Bourne-Again Shell (**bash**). However, if you are working on a Mac, your default might be the Z Shell (**zsh**). Any of these are are fine, but it is important to know which you are using for **unix profiles**.

If you want to know what Shell you're running you can just run `echo "$SHELL"`

Each time you start command-line or Terminal, your Shell will read specific files known as Shell initialization files to build your unix profile. The initialization files contain custom Shell configuration settings and can be found in your home directory (`cd ~`). The files are hidden (they begin with a period), but you can see hidden files with `ls -a`. 

If you are running **bash**, you will want to make your profile in **`.bash_profile`**. If you are running **zsh**, you will want to make your profile in **`.zshrc`**. If these files don't already exist, you may have to create them. You can use `nano` to both create and edit these files. 

``` bash
nano ~/.bash_profile
#OR
nano ~/.zshrc
```

There's lots of things you can put in your unix profile. For example, you can define shortcuts to do certain commands faster. Instead of typing `ls -lah`, why not just type `ll`. These are known as **aliases**. You can also define **functions** to do more complex things, like provide multiple arguments. An good explaination of alias and functions is given in the [Scripting OS X Blog](https://scriptingosx.com/2017/05/configuring-bash-with-aliases-and-functions/). 

Perhaps the most important thing in your unix profile is your **`$PATH`**, which is a variable containing a list of all possible directories where installed softwares may exist. Let's look at an example `.bash_profile` (this may be a good starting point for you):

``` bash
# Bash Profile

# User Profile
PS1="[Rhett@Macbook: \W] $ "

# SSH Profiles
alias remote="caffeinate ssh username@login.remote.server.edu"

# Alias Shortcuts
alias ll="ls -lah"
alias d="conda deactivate"
alias bio="conda activate bio"
alias envs="conda info --envs"

# Functions
mkcd ()
{
	mkdir -p -- "$1" && cd -P -- "$1"
}

# PATH
export PATH=/usr/local/bin:/usr/bin:/usr/sbin:/bin:/sbin
export PATH=$PATH:~/path/to/bin:~/path/to/bin/scripts
```

When I start a terminal/shell, **bash** will interpret my `.bash_profile`. It will tell it to print `[Rhett@Macbook: ~] $` at the beginning of each line. It will define each `alias` and `function` as a shortcut for specific commands. Finally, it will export my `$PATH` as: 
``` bash
/usr/local/bin:/usr/bin:/usr/sbin:/bin:/sbin:~/path/to/bin:~/path/to/bin/scripts
```
This string represents the location of several folders on my computer where scripts, commands, and softwares can be found. We haven't actually installed any softwares yet, but the next section will talk about a software management system that automatically adds softwares to your `$PATH`.

> **NOTE**
>
>`.zshrc` will look very similar - if not identical - to `.bash_profile`. If you are interested in moving from `bash` into the more feature-rich `zsh`, I recommend looking into:
>
>- [Scripting OS X: Moving to zsh](https://scriptingosx.com/2019/06/moving-to-zsh/)
>- [Oh My Zsh!](https://ohmyz.sh/)

## Software Management

[**Anaconda**](https://www.anaconda.com/distribution/#download-section) is a package/software management tool. I highly recommend you install Anaconda to make your life easier and avoid having to install many softwares by hand. 

When you install Anaconda, it will automatically add some script to your unix profile that will enable the Shell to find softwares installed by Anaconda. The reason Anaconda is so great is that software often do not play nicely with one another. For example, while one software may want to use `python v2.7`, another may want to use `python v3.6`. Anaconda allows you to create self-contained environments within which everything works nicely! If something needs a different version of a software…simply create a new environment! [Anaconda cheatsheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)

After you install Anaconda, we can create some environments for later use. Here’s how you do that:

``` bash
# The first step you only need to do once after installation.
# These lines set up what "channels" anaconda should in for software
# Configure conda to look in certain channels for packages
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels defaults

# Install a couple very useful packages in your base conda environment
conda install wget
conda install parallel

# Create a bioinformatics environment
conda create -n bio # Create an environment named "bio"
conda activate bio # Activate your "bio" environment
conda install -y biopython bamtools bedtools blast bowtie2 bwa cd-hit emboss fastqc gatk4 jellyfish parallel pear picard pigz rsem samtools sra-tools trim-galore
conda deactivate # Deactivate your environment to exit

# Create a Trinity environment
conda create -n trinity_env trinity parallel

# Create a ToxCodAn environment
conda create -n toxcodan_env python=3.7 biopython perl perl-bioperl perl-mce blast hmmer parallel

# Create a ChimeraKiller environment
conda create -n chimerakiller_env python=3.6 biopython bwa samtools bedtools picard pandas matplotlib scipy pysam gatk4 pathos parallel
```

Notice that I created a `bio` environment as well as some separate environments. This is because softwares often have different dependencies. For example, [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) requires a different version of [samtools](http://www.htslib.org/) than most other softwares installed in the `bio` environment. Putting softwares in their own environment ensures that they won’t interfere with each other. 

But how do you know when to create a new environment vs. just add to a pre-existing environment? Anaconda will generally warn you if it needs to change/upgrade/downgrade an existing software in that environment and will ask for confirmation to do this before continuing. When in doubt, create a new environment. You can create as many environments as you want and if something isn’t working properly, just create a new environment for it or remove the old environment and start over. All you have to do is remember to `conda activate` the correct environment before running an analysis.

I recommend to always try Googling "[`conda install {software}`](https://lmgtfy.app/?q=conda+install+trim+galore)" first; however, sometimes it is not possible to install things via Anaconda and you will have to manually download, install, and add the install location to your `$PATH`. I recommend creating a bin folder in an easily accessible location and installing all softwares in this folder. For example, I set up my bin folder in `~/Dropbox/bin` so that it can be accessed by any of my computers with Dropbox connected. You just need to make sure the `$PATH` is located in you unix profile on that computer.

Below, we download [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn), [CodAn](https://github.com/pedronachtigall/CodAn.git), and [signalp](https://services.healthtech.dtu.dk/software.php) and add them to our `$PATH`. I also recommend running the tutorial for each software to ensure it is installed properly.

``` bash
# Git clone the ToxCodAn repository and add to your PATH:
git clone https://github.com/pedronachtigall/ToxCodAn.git
echo "export PATH=\$PATH:$PWD/ToxCodAn/bin" >> ~/.bash_profile

# Git clone the CodAn repository and add to your PATH:
git clone https://github.com/pedronachtigall/CodAn.git
echo "export PATH=\$PATH:$PWD/CodAn/bin" >> ~/.bash_profile

# Download the SignalP-4.1, decompress and add it to your PATH:
tar -xvzf signalp-4.1g.Linux.tar.gz
echo "export PATH=\$PATH:$PWD/signalp-4.1/" >> ~/.bash_profile

# Source your .bash_profile to activate changes
source ~/.bash_profile
```
