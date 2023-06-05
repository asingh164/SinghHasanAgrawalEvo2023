

#Set up directory, import packages
setwd('/Users/asad/Drive/Desktop/EEB498/genome_scan')
library(tidyverse)
library(qvalue)

p_nogrm_true <- read_csv('p_nogrm_true.csv', col_names = FALSE)


    ### Get no. of blocks w/ p-values < 0.05, 0.01, 0.001 for each trait ###

#Choose trait of interest
trait = 'male.mean.relative.fitness.cage'

p_df <- p_nogrm_true %>%
  filter(p_nogrm_true[1] == trait)

temp <- as.numeric(p_df[1,3:1683])

length(temp[temp < 0.05]) # Get no. of blocks p < 0.05
length(temp[temp < 0.01]) # Get no. of blocks p < 0.01
length(temp[temp < 0.001]) # Get no. of blocks p < 0.001 

