# Substrate-based kinase activity inference


# load R packages
library(tidyverse)
library(parallel)


# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")


# load paths to input and output files
#args <- commandArgs(trailingOnly = TRUE)
args <- c("./data/colon_cancer_TvsN/ttop_tumourVsNormal.csv", "./output/files/colon_cancer_kinase_activities.txt")


# load input data (differential expression table)
diff_expr <- read_csv(file = args[1]) %>%
  separate(col = "ID", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(gene, position, residue, sample = t)


# load kinase-substrate list
# prepare kinase-substrate lists for inference
ks_lists <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct()


# infer kinase activities through a z-test without weights
kin_activities <- inferKA_W(
  ks_lists = ks_lists,
  phospho = diff_expr,
  with_weights = F,
  multi_core = F)

kin_activities <- kin_activities %>%
  select(-sample) %>%
  rename(substrates = n, activity = log10P)

write_tsv(kin_activities, args[2])
