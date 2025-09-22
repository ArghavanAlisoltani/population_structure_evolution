# install.packages(c("readr","readxl","dplyr","tibble"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicRanges","IRanges"))
setwd("~/Desktop/OSU_projects/conifers/LP/FSTcalc/provenances/100_greedy/fst_per_site")

# --- Inputs you change ---

library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(GenomicRanges)
library(data.table)

fstall <- data.frame(fread("fst_per_site_merged_rowmax_lt0.15.csv",sep=",") )
anno <- "~/Desktop/OSU_projects/conifers/LP/vcf_v1/richness/unimputed_all_renamed_1ab/SNP_nearest_annotation_ALL_snps.tsv"
P34<-data.frame(fread("PROC3_vs_PROC4.weir.fst"))
P35<-data.frame(fread("PROC3_vs_PROC5.weir.fst"))
P45<-data.frame(fread("PROC4_vs_PROC5.weir.fst"))

anno <- data.frame(fread(anno))
P34_top<-fstall[fstall$PROC3_vs_PROC4>=0.25,]
P34_top<-fstall[fstall$PROC3_vs_PROC4>=0.25,]

fstall$scaff_pos<-paste(fstall$CHROM,fstall$POS,sep = ":")
anno$scaff_pos<-paste(anno$X.CHROM,anno$POS,sep = ":")
merged<-merge(fstall,anno, by="scaff_pos", all.x=T)
fwrite(merged, "All_fst_merged_with_anno.tsv", sep = "\t")

