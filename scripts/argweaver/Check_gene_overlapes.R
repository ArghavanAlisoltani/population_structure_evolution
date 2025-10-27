#load libraries
library(data.table)
library(readxl)

#set directory and import data
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/all_tmrca/all_tmrca")
df<-data.frame(fread("fisher_mrna_out_v5/tmrca_with_mrna.tsv")) #start of the header
summary(df$mean_tmrca)
df_ordered<-df[order(df$mean_tmrca),]
df50<-df_ordered[c(1:562609),]
filt<-df[df$overlap_mrna>0,]
len2<-df[df$overlap_mrna_n>1 & df$seg_length<100,]
fwrite(filt, "all_overlaped.txt")

hist(df$mean_tmrca)
sum(df75$seg_length)
sum(df$seg_length)
sum(df50$seg_length)

nrow(df)*0.5


####################
#snp-wise
####################

snps <- read_excel("snps_annotated_tmrca_te_mrna.xlsx",
                             sheet = 2)
tbl<-data.frame(table(snps$tmrca_group))






