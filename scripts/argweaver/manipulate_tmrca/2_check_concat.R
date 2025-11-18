library(data.table)
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/oct_27_2025")
df1<-data.frame(fread("all.tmrca.txt",header = F,sep="\t"))
names(df1)<-c("CHROM","start_tmrca","end_tmrca",
              "mean_tmrca", "Low_CI", "UP_CI")
df1$seg_length<-(df1$end_tmrca-df1$start_tmrca)+1
df1no0<-df1[df1$mean_tmrca>0,]
df10<-df1[df1$mean_tmrca==0,]
summary(df1$mean_tmrca)
summary(df1$seg_length)
summary(df1no0$mean_tmrca)
summary(df1no0$seg_length)
sum(df1no0$seg_length)
sum(df1$seg_length)

df1test<-df1[c(140000:140001),]

tbl<-data.frame(table(df1$CHROM))
