options(stringsAsFactors = FALSE)

library("doBy")
library("RColorBrewer")
library("ggplot2")
library("pvclust")
library("pheatmap")
library("igraph")
library("plyr")
library("reshape")
library("picante")
library("pvclust")
library("gdata")

TraitData = read.csv("Trait_FullMatrix.csv", header = TRUE, check.names = F, row.names = 1)
IsolationData = read.csv("Isolation_FullMatrix.csv", header = TRUE, check.names = F, row.names = 1)
GeneraList = read.csv("GeneraList.csv", header = TRUE, check.names = F)

TraitData = TraitData[,-which(colnames(TraitData) == "Glucose")]

# Remove anything that is not a Ascomycetous fungus
TraitData$generaCol = rep(0, nrow(TraitData))
for(i in 1:nrow(TraitData)){
  TraitData[i,"generaCol"] = strsplit(rownames(TraitData)[i], " ")[[1]][1]
}
GeneraDrop = GeneraList[which(GeneraList[,2] == "B" | GeneraList[,2] == "NY"),1]
TraitData = TraitData[-which(TraitData[,"generaCol"] %in% GeneraDrop),]
TraitData = TraitData[-which(TraitData[,"generaCol"] == "Protomyces"),]
TraitData=TraitData[,-which(colnames(TraitData) == "generaCol")]

TreatNAs = data.frame(Treatment = character(), Number = numeric()) # Creates an empty dataframe to quantify the number of times a treatment was NA for the yeast
TREATMENT = 1 # Designates number for column to avoid string matching
NUMBER = 2 # Designates number for column to avoid string matching
for(i in 1:ncol(TraitData)){ # Start a loop to calculate the number NAs for a treatment
  TreatNAs[i, TREATMENT] = names(TraitData)[i] # Adds treatment name to row i of dataframe
  TreatNAs[i, NUMBER] = length(which(is.na(TraitData[,i]) == TRUE)) # Quantifies and adds the number of NAs to row i of the dataframe
} # End loop
AverageTests = mean(TreatNAs[,2])
KeepTreats = TreatNAs[which(TreatNAs[,"Number"] < 100),"Treatment"] # Determines which treatments have less than 100 NAs
TraitData = TraitData[,which(colnames(TraitData) %in% KeepTreats)] # Keeps treatments with less than 100 NAs

SpeciesNAs = data.frame(Species = character(), Number = numeric()) # Creates an empty dataframe to quantify the number of times a treatment was NA for the yeast
SPECIES = 1 # Designates number for column to avoid string matching
NUMBER = 2 # Designates number for column to avoid string matching
for(i in 1:nrow(TraitData)){ # Start a loop to calculate the number NAs for a treatment
  SpeciesNAs[i, SPECIES] = rownames(TraitData)[i] # Adds treatment name to row i of dataframe
  SpeciesNAs[i, NUMBER] = length(which(is.na(TraitData[i,]) == TRUE)) # Quantifies and adds the number of NAs to row i of the dataframe
} # End loop
AverageTests = mean(SpeciesNAs[,2])

rownames(IsolationData) = trimws(rownames(IsolationData), which = "right")

