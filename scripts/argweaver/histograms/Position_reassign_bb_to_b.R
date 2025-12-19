## Your data frame:
## mergedfinal
library(data.table)
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/NOV_27_2025/")
mergedfinal<-data.frame(fread("all.tmrca.txt",header = F,sep="\t"))
names(mergedfinal)<-c("CHROM","start_tmrca","end_tmrca",
              "mean_tmrca", "Low_CI", "UP_CI")
# Safety: ensure numeric (won’t hurt if already numeric)
mergedfinal$start_tmrca <- as.numeric(mergedfinal$start_tmrca)
mergedfinal$end_tmrca   <- as.numeric(mergedfinal$end_tmrca)

# Targets (exactly as provided)
targets <- c(
  "scaffold_6b","scaffold_7b","scaffold_8b","scaffold_9b","scaffold_10b",
  "scaffold_1ab","scaffold_1bb","scaffold_2b","scaffold_3b","scaffold_4b","scaffold_5b"
)

# Rows to modify
idx <- mergedfinal$CHROM %in% targets

# Add 1,000,000,000 to tmrca start/end on those rows
mergedfinal$start_tmrca[idx] <- mergedfinal$start_tmrca[idx] + 1e9
mergedfinal$end_tmrca[idx]   <- mergedfinal$end_tmrca[idx]   + 1e9

# Drop only the *last* 'b' from CHROM for those rows
mergedfinal$CHROM[idx] <- sub("b$", "", mergedfinal$CHROM[idx])

test<-mergedfinal[c(140000:140001),]

  
mergedfinal$seg_length<-(mergedfinal$end_tmrca-mergedfinal$start_tmrca)+1
nrow(mergedfinal)

fwrite(mergedfinal,"all_tmrca_corrected_position.tsv",sep="\t")


