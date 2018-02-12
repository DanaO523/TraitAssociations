options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("picante")
library("ape")
library("nnls")
library("phytools")
library("caper")
library("geiger")


#FASTA=read.fasta("D1D2_12-9-16.FASTA", seqtype = "DNA", as.string = TRUE)
#USENAMES=names(FASTA)
#SET ALL SEQUENCES TO CAPITAL IN NOTEPAD++
#FASTA=read.fasta("D1D2_12-9-16.FASTA", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
#names(FASTA) = USENAMES
#write.fasta(FASTA, names = USENAMES, as.string = TRUE, file.out = "AllSeq.fas")

TraitData = read.csv("traits_2018-01-10.csv", header = TRUE, row.names = 1)
Fix=read.csv("FixSpecies.csv", header = T)

tree<-read.tree(file="RAxML_bestTree_2017-01-11.result")
#dichotomousphylogeny <- di2multi(tree, random = FALSE, tol = 0.01)
dichotomousphylogeny<- multi2di(tree, random = FALSE, tol = 0.01)
rootedtree <- root(dichotomousphylogeny, "Filobasidiella_neoformans")
#rootedtree = dichotomousphylogeny
TipLabels = rootedtree$tip.label

for(i in 1:nrow(Fix)){
  ROWS =which(rownames(TraitData) == Fix[i,1])
  if(length(which(rownames(TraitData) == Fix[i,2])) == 0)
    rownames(TraitData)[ROWS] = Fix[i,2]
}


MissingData = rownames(TraitData)[-which(rownames(TraitData) %in% TipLabels)]
TraitData2=TraitData[-which(rownames(TraitData) %in% MissingData),]

SpeciesNames = TipLabels[-which(TipLabels %in% rownames(TraitData))]
rootedtree = drop.tip(rootedtree, SpeciesNames)

SpeciesOrder = rootedtree$tip.label
TraitOrdered_m = matrix(0, nrow = nrow(TraitData2), ncol = ncol(TraitData2))
colnames(TraitOrdered_m) = colnames(TraitData2)
TraitOrdered_m = data.frame(TraitOrdered_m, check.names = F)
for(i in 1:length(SpeciesOrder)){
  ROWS=which(rownames(TraitData2) == SpeciesOrder[i])
  TraitOrdered_m[i,] = TraitData2[ROWS,]
}
TraitOrdered_m$Species = SpeciesOrder
rownames(TraitOrdered_m) = SpeciesOrder

TraitOrdered_m2 = TraitOrdered_m[,-which(colnames(TraitOrdered_m) == "Species")]
TraitSums = data.frame(colSums(TraitOrdered_m2))
TraitSums$Traits = rownames(TraitSums)
colnames(TraitSums) = c("TraitCount", "Traits")

TC= ggplot(TraitSums, aes(x = Traits, y = TraitCount))+geom_point()+geom_hline(yintercept = 0, colour = "gray60", lty = 2)
TC = TC + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.2))
ggsave(paste("TraitCount_", Sys.Date(), ".pdf",sep = ""),useDingbats = F)


#ProbTraits = rownames(TraitSums)[which(TraitSums[,1] == 0 | TraitSums[,1] == 561)]
#TraitOrdered_m3 = TraitOrdered_m2[,-which(colnames(TraitOrdered_m2) %in% ProbTraits)]
TraitOrdered_m3 = TraitOrdered_m2

SPECIES = rownames(TraitOrdered_m3)
TRAITS = colnames(TraitOrdered_m3)

Trait_df = data.frame(Traits = character(), Species = character(), PresAbs = numeric())
k = 1
for(i in 1:length(TRAITS)){
  for(j in 1:length(SPECIES)){
    Trait_df[k,1] = colnames(TraitOrdered_m3)[[i]]
    Trait_df[k,2] = rownames(TraitOrdered_m3)[[j]]
    Trait_df[k,3] = TraitOrdered_m3[j,i]
    k = k+1
  }
}

# PhyloD_list = list()
# for(i in 1:length(TRAITS)){
#   TEMPTRAIT = Trait_df[which(Trait_df$Traits == TRAITS[[i]]),]
#   CompData = comparative.data(data = TEMPTRAIT, phy = rootedtree, names.col = `Species`, vcv = TRUE)
#   TEMP =phylo.d.subset(data = CompData, phy = rootedtree, binvar = PresAbs)
#   PhyloD_list[[i]] = TEMP
# }

PhyloD_list2 = list()
for(i in 1:length(TRAITS)){
  TEMPTRAIT = Trait_df[which(Trait_df$Traits == TRAITS[[i]]),]
  CompData = comparative.data(data = TEMPTRAIT, phy = tree, names.col = `Species`, vcv = TRUE)
  TEMP =phylo.d(data = CompData, phy = tree, binvar = PresAbs)
  PhyloD_list2[[i]] = TEMP
}
names(PhyloD_list2) = TRAITS

PhyloD_df = data.frame(Trait = character(), D = numeric(), P_Random = numeric(), P_Brownian = numeric())

for(i in 1:length(PhyloD_list2)){
  PhyloD_df[i,1] = names(PhyloD_list2)[i]
  PhyloD_df[i,2] =  PhyloD_list2[[i]][1]
  PhyloD_df[i,3] =  PhyloD_list2[[i]][2]
  PhyloD_df[i,4] =  PhyloD_list2[[i]][3]
}

TraitPresent = colSums(TraitOrdered_m3)
TraitPresent_df = data.frame(Trait = names(TraitPresent), Count = TraitPresent)
PhyloD_df = merge(PhyloD_df, TraitPresent_df)
write.csv(PhyloD_df, file = "PhyloD_df2.csv")

