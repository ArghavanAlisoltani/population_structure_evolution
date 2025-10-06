library(vroom)
setwd("~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly")
my_data <- data.frame(vroom("~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/readable.TEanno.cds.scaf1ab.tsv",
                 col_names = F))
class(my_data)
names(my_data)<-c("seqid", "source", "sequence_ontology", "start",
             "end", "score", "strand", "phase", "attributes")

head(my_data[c(2000000:2000010),c(1:8)])
Scaff_241ab<-c(
  "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5",
  "scaffold_6", "scaffold_7", "scaffold_8", "scaffold_9",
  "scaffold_10", "scaffold_11", "scaffold_12", "scaffold_13",
  "scaffold_14", "scaffold_15", "scaffold_16", "scaffold_18",
  "scaffold_19", "scaffold_20", "scaffold_27", "scaffold_31",
  "scaffold_44", "scaffold_186", "scaffold_528", "scaffold_1a",
  "scaffold_1b"
)

my_data<-my_data[my_data$seqid %in% Scaff_241ab,-c(8:9)]
fwrite(my_data,"~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv",
         sep="\t")
         
TE<-data.frame(fread("~/Desktop/OSU_projects/conifers/LP/lodgepole_pine_assembly/col8_readable.TEanno.cds.scaf1ab.tsv",
                     sep="\t", header=T
                     #skip = 4500000, 
                     #nrow=4700000
                     ))
table<-data.frame(table(TE$sequence_ontology))
fwrite(table,"sequence_ontology.tsv",sep = "\t")



         