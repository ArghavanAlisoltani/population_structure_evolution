library(data.table)
setwd("~/Desktop/OSU_projects/conifers/LP/ARGweaver/all_tmrca/all_tmrca/py_tmrca")
df<-data.frame(fread("All_tmrca_60_run.txt", sep="\t", header=F))
summary(df$V4)
df_non_zero<-df[df$V4>0,]
df_zero<-df[df$V4==0,]
summary(df_non_zero$V4)
boxplot(df_non_zero$V4~df_non_zero$V1, outline=F)

# Create the violin plot
# Make sure you have the ggplot2 library installed and loaded
# install.packages("ggplot2")
library(ggplot2)

# Create the violin plot
ggplot(df_non_zero, aes(x = V1, y = V4)) +
  #ylim(c(0,20000))+
  
  # Add the violin layer with a single fill color
  geom_violin(trim = FALSE, fill = "gray") +
  
  # Add a narrow box plot and hide outliers
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  
  # Add a point for the mean value
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "black") +
  
  # Add labels and a title for clarity
  labs(
    title = "TMRCA Distribution by Group",
    x = "Scaffold",
    y = "TMRCA"
  ) +
  
  # Use a clean theme
  theme_classic() +
  
  # Rotate x-axis labels to 90 degrees
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
