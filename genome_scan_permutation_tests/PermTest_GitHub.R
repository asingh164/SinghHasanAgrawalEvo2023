
#Goals:
# - Run permutation tests on 8 traits (SA-vial, SC-vial, SA-cage, SC-cage, female-vial, male-vial, female-cage, male-cage)
#   and collect:
#       (i) p-value and effect size of relatedness matrix
#       (ii) # of blocks with p < 0.001, p < 0.01, p < 0.05, p < 0.01, Q < 0.1, Q < 0.2
#       (iii) For each of the significance thresholds, the number of clusters (moderate and conservative)

#Packages
library(dplyr)
library(tidyr)
library(readr)

#Load in Data
setwd('/Users/asad/drive/Desktop/EEB498/AsadProject2020-main')

pheno.data <- read.csv('data_summary_vial_cages.csv')
HMM.matrix <- read.delim('HMM_matrix.txt')


######################################################################################
###   In house Genome Scan - Permutation Tests      ###
######################################################################################

#loading packages
require(coxme)
require(psych)
require(qvalue)
require(gdata)

#Loading functions


######Prepping files for lmekin gwas ######

args <- commandArgs(trailingOnly = TRUE)

#Prepping output files to collect block-by-block grm model p-values
#Preparing Phenotype Data
pheno.data = pheno.data[order(pheno.data$DSPR.stock.number),] #Put RILIDs in numeric order


list.of.positions = unique(HMM.matrix$row.id)[seq(1, length(unique(HMM.matrix$row.id)), 7)] #List of unique block IDs in the genome upon which to conduct the genome scan
threshold = 0.1  #This is the threshold we will use to exclude rils below


#List of traits that I want to perform the genome scan with (currently the 8 traits)
list.of.traits = (names(pheno.data)[2:9])

for (a in 1:1000){ 
  
  message(paste('starting permutation number', a)) 
  
  #Fake phenotype data
  fake_pheno <- pheno.data[2:9] %>%
    mutate(fake.id = sample(pheno.data$DSPR.stock.number, replace = FALSE)) %>%
    select(fake.id, everything()) %>%
    arrange(fake.id)
  colnames(fake_pheno)[1] <- 'DSPR.stock.number'

  
  #Setting up the output dataframes
  p.nogrm.scan = as.data.frame(matrix(NA, ncol = length(list.of.positions) + 1, nrow = length(list.of.traits)))
  p.nogrm.scan$V1 <- list.of.traits
  colnames(p.nogrm.scan) <- c('traits', as.character(list.of.positions))
  
  d.nogrm.scan = as.data.frame(matrix(NA, ncol = length(list.of.positions) + 1, nrow = length(list.of.traits)))
  d.nogrm.scan$V1 <- list.of.traits
  colnames(d.nogrm.scan) <- c('traits', as.character(list.of.positions))
  
  
  
  
  
  trait.number = 1
  
  for (t in list.of.traits){
    
    t = list.of.traits[1]
    output.col = 0

    #Subsetting the HMM matrix for RILs that I assayed for that specific trait
    pheno.data.subset = subset(fake_pheno, fake_pheno[,t] != "NA")
    HMM.matrix.subset.rils = HMM.matrix[HMM.matrix$RILID %in% pheno.data.subset$DSPR.stock.number, ]
    
    #Genome scan starts here
    for (i in 1:length(list.of.positions)){
      
      i = 1
      output.col = output.col + 1
      HMM.matrix.subset = subset(HMM.matrix.subset.rils, HMM.matrix.subset.rils$row.id == list.of.positions[i])
      geno.prob.sums = apply(HMM.matrix.subset[, 4:11], 2, function(x) sum(x))  #Determine the cumulative probability of each genotype at the specific block
      genotypes.to.exclude = subset(geno.prob.sums, geno.prob.sums <= 5)  #This is a subset of the genotypes that reach a cumulative probability of less than 5% across all RILs
      
      #Now I want to remove RILs in which the probability that a RIL is of the genotypes to be removed is greater than the threshold probability
      HMM.matrix.subset.for.lm = HMM.matrix.subset
      for (j in names(genotypes.to.exclude)){
        HMM.matrix.subset.for.lm = HMM.matrix.subset[!(HMM.matrix.subset[,j] > threshold),]
      }
      merged.df = merge(HMM.matrix.subset.for.lm, fake_pheno, by.x = "RILID" , by.y = "DSPR.stock.number")
      
      ##Now that we have our HMM list, we want to update our model statement
      genotypes.to.include = names(subset(geno.prob.sums,  geno.prob.sums > 5)) #This is a subset of the genotypes that reach a cumulative probability of more than 5% across all RILs
      number.of.factors = length(genotypes.to.include) - 1
      
      #Determining the elements that will be part of the model statement depending on which genotypes will be included in the model
      formula.vector = c()
      for (k in 1:(length(genotypes.to.include) - 1)){
        formula.vector[k] = paste("merged.df[,'", genotypes.to.include[k], "'] ", sep = "")
      }
      
      #Now  will dynamically update the model statement based on the genotypes included as factors in the model
      response.variable = paste(as.character(t), "~", sep = "")
      model.statement = as.formula(paste(response.variable, paste(formula.vector, collapse = "+"), "+ (1|RILID)"))
      model.statement.null = as.formula(paste(response.variable, "1 + (1|RILID)"))
      
      #Fitting LM without GRM and extracting summary statistics
      model.ml = lm(model.statement, data=merged.df)
      model.ml.null = lm(model.statement.null, data=merged.df)
  
      
      #Extracting p-value and test statistic for anova between null and alternate model w/o GRM as random effect 
      p.nogrm.scan[trait.number, output.col + 1] = anova(model.ml.null, model.ml)[2,6]
      d.nogrm.scan[trait.number, output.col + 1] = anova(model.ml.null, model.ml)[2,5]
      
      message(paste('finished block number '), output.col)
      
      
    }
    
    message(paste('finished '), list.of.traits[trait.number])
    trait.number = trait.number + 1
    
    
    
    
  }
  
  #Define permutation number in output files
  p.nogrm.scan$permutation <- a
  d.nogrm.scan$permutation <- a
  
  #Arrange output data for output file
  p.nogrm.scan <- select(p.nogrm.scan, permutation, traits, everything())
  d.nogrm.scan <- select(d.nogrm.scan, permutation, traits, everything())
  
  #Create output files and append subsequent permutation data to output file
  write_csv(p.nogrm.scan, path = paste('nogrm_p', as.character(args), '.csv', sep = ''), quote = FALSE, append = TRUE,
            col_names = FALSE)
  
  write_csv(d.nogrm.scan, path = paste('nogrm_d', as.character(args), '.csv', sep = ''), quote = FALSE, append = TRUE,
            col_names = FALSE)
  
}






