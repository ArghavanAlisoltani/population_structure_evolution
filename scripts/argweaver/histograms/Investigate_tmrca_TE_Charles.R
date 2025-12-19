setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025")
library(data.table)
TE<-fread("~/Desktop/Shiny_Apps/BFF_AFIRMS/col8_readable.TEanno.cds.scaf1ab.tsv")
unique<-data.frame(table(TE$sequence_ontology))
TE1<-TE[TE$seqid=="scaffold_1b" & TE$start==447420149 | TE$start==447420154,]
TE1$length<-(TE1$end-TE1$start)+1
duplicated<-data.frame(table(TE$start,TE$sequence_ontology))
tmrca<-fread("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/all_tmrca_corrected_position.tsv")
tmrca1<-tmrca[tmrca$CHROM=="scaffold_1b"& tmrca$start_tmrca>=447420110 & tmrca$end_tmrca<=447437260]
mean(tmrca1$mean_tmrca)
median(tmrca1$mean_tmrca)

TE_tmrca<-fread("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/tmrca_fast_bins_full_v2_mybins_v5/te_tmrca_all.tsv")
TE_tmrca1<-TE_tmrca[TE_tmrca$seqid=="scaffold_1b"& TE_tmrca$start>=447420110 & TE_tmrca$end<=447437260]

fwrite(tmrca1,"tmrca_data_for_Charles/TMRCA_scaffold1b_447420110_447437260.txt", sep="\t")

