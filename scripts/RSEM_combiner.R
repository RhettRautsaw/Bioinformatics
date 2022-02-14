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
# Combine RSEM Results for multiple individuals
# 
#=============#
##### USE #####
#=============#
#
# RSEM_combiner.R -p .genes.results -o species_expression.txt -g TOXIN -c 0.95
#
#======================#
##### REQUIREMENTS #####
#======================#

packages = c("argparse","readr", "stringr", "tidyr", "dplyr", "data.table")

# Install packages not yet installed
installed_packages = packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  print("Installing required packages")
  install.packages(packages[!installed_packages], repos="https://cloud.r-project.org")
}


#=========================#
##### SETUP ARGUMENTS #####
#=========================#

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# create parser object
parser = ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option
parser$add_argument("-p", default=".genes.results", 
                    help="String pattern to search for to find relevent files [default: \"%(default)s\"]")
parser$add_argument("-o", default="species_expression.txt", 
                    help="Name of output file [default: \"%(default)s\"]")
parser$add_argument("-g", default="NA", 
                    help="Optional: String pattern to search for in gene names to create group column [default: \"%(default)s\"]")
parser$add_argument("-c", default=0, 
                    help="Optional: Quantile cutoff. Group (-g) below this limit are flagged and placed in a list for you to remove from the transcriptome [default: \"%(default)s\"]")

#========================#
##### READ ARGUMENTS #####
#========================#

args = parser$parse_args()
# args$g="TOXIN"
# args$c=0.95

## MAKE LIST OF FILES
files = list.files(path=".", pattern=args$p)
names = gsub(args$p, "", files)

## READ FILE 1
merged = readr::read_delim(files[1], show_col_types = FALSE) %>% select(-`transcript_id(s)`, -FPKM)
setnames(merged, old = c("expected_count","TPM"), new = c(paste0("COUNT_",names[1]), paste0("TPM_",names[1]))) #

## MERGE ALL OTHER FILES
if(length(files)>1){
  for(i in 2:length(files)){
    tmp = readr::read_delim(files[i], show_col_types = FALSE) %>% select(-`transcript_id(s)`, -FPKM, -length, -effective_length)
    setnames(tmp, old = c("expected_count","TPM"), new = c(paste0("COUNT_",names[i]), paste0("TPM_",names[i]))) # 
    merged = merge(merged,tmp, by="gene_id")
  }
}

## SORT COLUMNS
merged = merged[ , c(1,2,3,order(names(merged)[c(-1,-2,-3)])+3)]

## CALCULATE AVERAGE TPM
TPM_columns = which(str_detect(names(merged),"TPM"))
merged$TPM_average = rowMeans(merged[,TPM_columns])
merged = merged %>% select(gene_id,length,effective_length,TPM_average,everything())

if(args$g!="NA"){
  ## Create Group Column & Sort
  merged$group = str_extract(merged$gene_id,args$g)
  merged$group = replace_na(merged$group,"OTHER")
  
  ## Create Cutoff Column & Sort
  cutoff = merged %>% group_by(group) %>% summarize(cutoff=quantile(TPM_average,as.numeric(args$c))) %>% 
    filter(group=="OTHER") %>% select(cutoff) %>% as.numeric()
  merged = merged %>% mutate(above_cutoff = ifelse(group==args$g, TPM_average>cutoff, NA)) %>%
    select(gene_id, group, length, effective_length, above_cutoff, TPM_average, everything()) %>%
    arrange(desc(group),gene_id)
  low_list = merged %>% filter(above_cutoff==FALSE) %>% select(gene_id)
  write_delim(low_list, "v0_low_expression.list", delim="\t", col_names = F, quote="none")
}

write_delim(merged, args$o, delim="\t", quote="none")




