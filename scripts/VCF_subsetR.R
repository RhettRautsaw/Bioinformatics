#!/usr/bin/env Rscript

# Copyright 2022 Rhett M. Rautsaw
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#=====================#
##### DESCRIPTION #####
#=====================#
#
# VCF_subsetR is designed to take a vcf/bcf file and subset it according to provided metadata
#
#=============#
##### USE #####
#=============#
#
# VCF_subsetR -i input.vcf.gz -m metadata.txt -f "Tissue == 'Host' & Year %in% 2006:2009 & Site =='WPP'"
#
#======================#
##### REQUIREMENTS #####
#======================#

cat("\n\nStarting VCF_subsetR\n\n\n")
cat("Checking for proper installations\n")

packages = c("argparse","readr","dplyr","stringr","rlang","parallel")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  print("Installing required packages")
  install.packages(packages[!installed_packages], repos="https://cloud.r-project.org")
}

invisible(system(paste0("vcftools --version")))
invisible(system(paste0("bcftools --version")))

suppressPackageStartupMessages(library("argparse"))
cpu = parallel::detectCores()

#=========================#
######## ARGUMENTS ########
#=========================#
parser = ArgumentParser(description="VCF_subsetR will subset a vcf file given filter criteria and associated metadata")

parser$add_argument("-i","--input",
                    type="character",
                    default="input.vcf.gz",
                    help="VCF/BCF(.gz) file to filter. Default: None")
parser$add_argument("-m","--metadata",
                    type="character",
                    default="metadata.txt",
                    help="Tab-delimited metadata with first column equal to sample names in the input vcf/bcf. Default: None")
parser$add_argument("-f","--filter",
                    type="character",
                    default="Tissue == \'Host\' & Year %in% 2006:2009 & Site ==\'WPP\'",
                    help="String for filtering database in dplyr format. Note you will need to escape quotation marks with a backslash. For example, \"Tissue == \\'Host\\' & Year %%in%% 2006:2009 & Site ==\\'WPP\\'\". Default: None")
parser$add_argument("-o","--out",
                    type="character",
                    default="VCF_subset",
                    help="Folder/Prefix in which to output results. Default: VCF_subset")
parser$add_argument("-c","--cpu",
                    type="integer",
                    default=cpu,
                    help="Number of processors to be used in each step. (default: %(default)s)")

args = parser$parse_args()

#=======================#
######### SETUP #########
#=======================#

input_path = normalizePath(args$i)
if(grepl(".vcf.gz", input_path)){
  type="--gzvcf"
}else if(grepl(".bcf|.bcf.gz", input_path)){
  type="--bcf"
} else {
  type="--vcf"
}

metad_path = normalizePath(args$m)
filtr = args$f
out_name = args$o
out_path = paste0(getwd(),"/",out_name)
out_files = paste0(out_path,"/",out_name)

cpus = args$cpu

cat("\n\nVCF_subsetR Parameters:\n")
cat(paste0("\t Input:\t\t", args$i,"\n"))
cat(paste0("\t Metadata:\t", args$m,"\n"))
cat(paste0("\t Filter:\t", args$f,"\n"))
cat(paste0("\t Output:\t", args$o,"\n"))
cat(paste0("\t Processors:\t", args$c,"\n"))

#=======================#
##### LOAD PACKAGES #####
#=======================#

cat(paste0("\n", Sys.time(), " ::: Loading packages and preparing directory :::\n"))
suppressWarnings(invisible(suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))))

#=======================#
######### CODE ##########
#=======================#

# Make output folder
if(!dir.exists(out_path)){
  dir.create(out_path)
}

# Read/Filter Metadata
cat(paste0(Sys.time(), " ::: Reading and Filtering Metadata :::\n"))
metad = read_delim(metad_path, delim="\t", show_col_types = FALSE)
metad_filtr = metad %>% filter(!!parse_expr(filtr))
metad_filtr_list = metad_filtr %>% pull(1)
if(length(metad_filtr_list)==0){
  quit()
}
write_delim(x = metad_filtr, 
            file = paste0(out_path,"/metadata.txt"),
            delim="\t")
write_lines(x = metad_filtr_list, 
            file = paste0(out_path,"/samples.list"))

# Subset VCF file to matching samples & BGZIP
cat(paste0(Sys.time(), " ::: Subsetting VCF/BCF :::\n"))
system(paste0("vcftools ", type, " ", input_path, " --keep ", out_path, "/samples.list --recode --recode-INFO-all --out ", out_files))
system(paste0("mv ", out_files, ".recode.vcf ", out_files, ".vcf"))
cat(paste0(Sys.time(), " ::: Zipping VCF :::\n"))
system(paste0("bgzip -@ ", cpus, " ", out_files,".vcf"))
cat(paste0(Sys.time(), " ::: Indexing VCF :::\n"))
system(paste0("bcftools index --threads ", cpus, " ", out_files,".vcf.gz"))

# Print Helpful Message
cat(paste0("

    ALL FINISHED! You now have a subset fasta! 

You may now want to do some additional filtering with vcftools.
I'd recommend the following settings to start with to filter for 
biallelic SNPs with 5-100x coverage and 95% completeness (or only 5% missing data):

  vcftools --gzvcf input.vcf.gz --recode --recode-INFO-all
      --min-alleles 2 --max-alleles 2 --minDP 5 --maxDP 100 --max-missing 0.95 
      --out output_prefix

Some other options to consider may include: --maf (0.01-0.05) --mac --hwe --min-meanDP --max-meanDP --minQ
You may also want to consider an iterative approach to SNP filtering as suggested by O'Leary et al. 2018
  See Table 2: https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14792
\n\n"))
