
TRAITS = colnames(traits)
PERMUTEDVAL = 10000

#traits = read.csv("traits.csv", header = TRUE)
#load("PermutedAssoc_list.rda")

############# This is to generate permuted data - Permuted Data gets loaded above ##############
PermutedData= list()
for(iter in 1:PERMUTEDVAL){
  PermutedData[[iter]] = list()
  tempMatrix <- randomizeMatrix(traits, null.model="independentswap", iterations=1000)
  tempMatrix = data.frame(tempMatrix, check.names = F)
  tempMatrix$Species = rownames(tempMatrix)
  PermutedData[[iter]]$Permutation <- iter
  PermutedData[[iter]]$Permutation = iter
  PermutedData[[iter]]$CompiledData <- tempMatrix
  print(iter)
}

save(PermutedData, file = "PermutedData.rda")
#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)] ###########
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)       ###########
#############################################################################################################

# Calculates trait/trait counts for permuted data #
PermutedAssoc_list = list()
for(perm in 1:PERMUTEDVAL){
  PermutedAssoc_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData
  tempDF=tempDF[,which(colnames(tempDF) %in% TRAITS)]
  tempMatrix_Pos = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_Pos) = colnames(tempDF)
  rownames(tempMatrix_Pos) = colnames(tempDF)
  tempMatrix_Neg = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
  colnames(tempMatrix_Neg) = colnames(tempDF)
  rownames(tempMatrix_Neg) = colnames(tempDF)
  for(i in 1:ncol(tempDF)){
    for(j in 1:ncol(tempDF)){
      POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1))
      NEG1 = length(which(tempDF[,i] == 1 & tempDF[,j] == 0 ))
      NEG2 = length(which(tempDF[,i] == 0 & tempDF[,j] == 1))
      tempMatrix_Pos[i,j] = POS
      tempMatrix_Neg[i,j] = NEG1+NEG2
    }
  }
  PermutedAssoc_list[[perm]]$Positive = tempMatrix_Pos
  PermutedAssoc_list[[perm]]$Negative = tempMatrix_Neg
  print(perm)
}
# Calculates trait/trait counts for observed data #
TraitCountPos_m = matrix(0, nrow = length(TRAITS), ncol = length(TRAITS))
rownames(TraitCountPos_m) = TRAITS
colnames(TraitCountPos_m) = TRAITS
for(i in 1:ncol(traits)){
  for(j in 1:ncol(traits)){
    POS=length(which(traits[,i] == 1 & traits[,j] == 1))
    TraitCountPos_m[which(rownames(TraitCountPos_m) == colnames(traits)[i]),which(colnames(TraitCountPos_m)== colnames(traits)[j])] = POS
  }
}

Assoc_Pos_df = data.frame(Trait_A = character(),
                          Trait_B = character(),
                          Observed = numeric(),
                          Expected = numeric(),
                          Percent = numeric(),
                          Count = numeric(),
                          Pvalue = numeric())
k = 1
for(i in 1:length(TRAITS)){
  for(j in 1:length(TRAITS)){
    ROWI = which(colnames(TraitCountPos_m) == TRAITS[[i]])
    COLJ= which(rownames(TraitCountPos_m) == TRAITS[[j]])
    
    tempObs = TraitCountPos_m[i,j]
    tempList = list()
    for(perm in 1:PERMUTEDVAL){
      tempDF1 = PermutedAssoc_list[[perm]]$Positive
      OBS = tempDF1[ROWI,COLJ]
      tempList[perm] = OBS
    }
    temp = unlist(tempList)
    N = sum(temp >= tempObs)
    p = zapsmall(binconf(N, length(temp), method = "exact"))
    
    Assoc_Pos_df[k, 1] = TRAITS[[i]]
    Assoc_Pos_df[k, 2] = TRAITS[[j]]
    Assoc_Pos_df[k, 3] = tempObs
    Assoc_Pos_df[k, 4] = mean(temp, na.rm = T)
    Assoc_Pos_df[k, 5] = ecdf(temp)(tempObs)
    Assoc_Pos_df[k, 6] = N
    Assoc_Pos_df[k, 7] = p[1]
    k = k + 1
  }
}


for(i in 1:nrow(Assoc_Pos_df)){
  TraitA = Assoc_Pos_df[i,"Trait_A"]
  TraitB = Assoc_Pos_df[i,"Trait_B"]
  ROWDROP = which(Assoc_Pos_df$Trait_A == TraitB & Assoc_Pos_df$Trait_B == TraitA)
  if(length(ROWDROP > 0)){
    Assoc_Pos_df[ROWDROP,] = 0
	}
}

Assoc_Pos_df = Assoc_Pos_df[-which(Assoc_Pos_df$Trait_A == 0),] 

Assoc_Pos_df$Sig = Assoc_Pos_df$Pvalue <= 0.05
Assoc_Pos_df$Correction = p.adjust(Assoc_Pos_df$Pvalue, method = "BH")

###### Determine negative associations between Trait pairs ######
TraitCountNeg_m = matrix(0, nrow = length(TRAITS), ncol = length(TRAITS))
rownames(TraitCountNeg_m) = TRAITS
colnames(TraitCountNeg_m) = TRAITS
for(i in 1:ncol(traits)){
  for(j in 1:ncol(traits)){
    NEG1 = length(which(traits[,i] == 1 & traits[,j] == 0))
    NEG2 = length(which(traits[,i] == 0 & traits[,j] == 1))
    TraitCountNeg_m[which(rownames(TraitCountNeg_m) == colnames(traits)[i]),which(colnames(TraitCountNeg_m)== colnames(traits)[j])] = NEG1+NEG2
  }
}


Assoc_Neg_df = data.frame(Trait_A = character(),
                              Trait_B = character(),
                              Observed = numeric(),
                              Expected = numeric(),
                              Percent = numeric(),
                              Count = numeric(),
                              Pvalue = numeric())
k = 1
for(i in 1:length(TRAITS)){
  for(j in 1:length(TRAITS)){
    ROWI = which(colnames(TraitCountNeg_m) == TRAITS[[i]])
    COLJ= which(rownames(TraitCountNeg_m) == TRAITS[[j]])
    
    tempObs = TraitCountNeg_m[i,j]
    tempList = list()
    for(perm in 1:PERMUTEDVAL){
      tempDF1 = PermutedAssoc_list[[perm]]$Negative
      OBS = tempDF1[ROWI,COLJ]
      tempList[perm] = OBS
    }
    temp = unlist(tempList)
    N = sum(temp >= tempObs)
    p = zapsmall(binconf(N, length(temp), method = "exact"))
    Assoc_Neg_df[k,1] = TRAITS[[i]]
    Assoc_Neg_df[k,2] = TRAITS[[j]]
    Assoc_Neg_df[k,3] = tempObs
    Assoc_Neg_df[k,4] = mean(temp, na.rm = T)
    Assoc_Neg_df[k,5] = ecdf(temp)(tempObs)
    Assoc_Neg_df[k,6] = N
    Assoc_Neg_df[k,7] = p[1]
    k = k + 1
  }
}

for(i in 1:nrow(Assoc_Neg_df)){
  TraitA = Assoc_Neg_df[i,"Trait_A"]
  TraitB = Assoc_Neg_df[i,"Trait_B"]
  ROWDROP = which(Assoc_Neg_df$Trait_A == TraitB & Assoc_Neg_df$Trait_B == TraitA)
  if(length(ROWDROP > 0)){
    Assoc_Neg_df[ROWDROP,] = 0
  }
}

Assoc_Neg_df = Assoc_Neg_df[-which(Assoc_Neg_df$Trait_A == 0),] 


Assoc_Neg_df$Sig = Assoc_Neg_df$Pvalue <= 0.05
Assoc_Neg_df$Correction = p.adjust(Assoc_Neg_df$Pvalue, method = "BH")

Assoc_Pos_df$Direction = rep("Random", nrow(Assoc_Pos_df))
Assoc_Neg_df$Direction = rep("Random", nrow(Assoc_Neg_df))

for(i in 1:nrow(Assoc_Pos_df)){
  if(Assoc_Pos_df[i,"Correction"] < 0.05)
    Assoc_Pos_df[i,"Direction"] = "Positive"
}

for(i in 1:nrow(Assoc_Neg_df)){
  if(Assoc_Neg_df[i,"Correction"] < 0.05)
    Assoc_Neg_df[i,"Direction"] = "Negative"
}

PooledAssc_df = data.frame(rbind(Assoc_Pos_df, Assoc_Neg_df))

AllAsscData = matrix(0, nrow = 1, ncol = ncol(PooledAssc_df))
colnames(AllAsscData) = colnames(PooledAssc_df)
AllAsscData = data.frame(AllAsscData)

k = 1

for(i in 1:length(TRAITS)){
  for(j in 1:length(TRAITS)){
    tempData = PooledAssc_df[which(PooledAssc_df$Trait_A == TRAITS[i] & PooledAssc_df$Trait_B == TRAITS[j]),]
    if(length(unique(tempData$Direction)) == 1){
      AllAsscData[k,] = tempData[1,]
      k = k + 1
    }else if(length(unique(tempData$Direction)) >1 & length(which(tempData$Direction == "Random")) == 1){
      AllAsscData[k,] = tempData[which(tempData$Direction != "Random"),]
      k = k + 1
    }else if(length(unique(tempData$Direction)) >1 & length(which(tempData$Direction == "Random")) > 1){
      AllAsscData[k, ] = tempData[1,]
      AllAsscData[k,"Direction"] = "Shit"
      k = k + 1
      AllAsscData[k, ] = tempData[2,]
      AllAsscData[k,"Direction"] = "Shit"
      k = k+1
    }
  }
}


#AllAsscData[which(AllAsscData$Direction == "Shit"),"Direction"] = "Random"
AllAsscData$TraitDiff = AllAsscData$Observed - AllAsscData$Expected
AllAsscData$AbsDiff = abs(AllAsscData$TraitDiff)
write.csv(AllAsscData, file = "PooledAssociation_df.csv")

PooledAssoc_df = AllAsscData 
k = 1
for(i in 1:nrow(AllAsscData)){
  PooledAssoc_df[k,] = AllAsscData[i,]
  k = k+1
  tempDF = AllAsscData[i,]
  tempDF[1,"Trait_A"] = AllAsscData[i,"Trait_B"]
  tempDF[1,"Trait_B"] = AllAsscData[i,"Trait_A"]
  PooledAssoc_df[k,] = tempDF
  k = k+1
}

NEGROW = which(PooledAssoc_df$Direction == "Negative")
RANDROW = which(PooledAssoc_df$Direction == "Random")

test  = PooledAssoc_df
for(i in 1:length(NEGROW)){
  test[NEGROW[[i]],"AbsDiff"] = (test[NEGROW[[i]],"AbsDiff"])*-1
}

for(i in 1:length(RANDROW)){
  test[RANDROW[[i]],"AbsDiff"] = (test[RANDROW[[i]],"AbsDiff"])*0
}

Assoc_m = matrix(0, nrow = length(TRAITS), ncol = length(TRAITS))
colnames(Assoc_m) = TRAITS
rownames(Assoc_m) = TRAITS

for(i in 1:nrow(test)){
  Assoc_m[which(rownames(Assoc_m) == test[i,"Trait_A"]),which(colnames(Assoc_m) == test[i,"Trait_B"])] = test[i,"AbsDiff"]
}

fit = pvclust(Assoc_m, method.dist = "euclidean", method.hclust = "ward.D", nboot = 1000, store = F)
clusters= pvpick(fit, alpha = 0.95,max.only = F)

pdf("TxT_dendro.pdf")
plot(fit)
pvrect(fit, alpha = 0.95, max.only = TRUE)
dev.off()

ord <- hclust( dist(Assoc_m, method = "euclidean"), method = "ward.D" )$order

pd <- as.data.frame( Assoc_m )
pd$Trait_A <- rownames(Assoc_m)
pd.m <- melt( pd, id.vars = "Trait_A", variable.name = "Trait_B" )
colnames(pd.m) = c("Trait_A", "Trait_B", "Diff")

pd.m$Trait_B <- factor( pd.m$Trait_B, levels = colnames(Assoc_m)[ord])
pd.m$Trait_A <- factor( pd.m$Trait_A, levels = rownames(Assoc_m)[ord])

# Melt the correlation matrix
colnames(pd.m) = c("Trait_A", "Trait_B", "Difference")
Assoc_df = pd.m
Assoc_df$Direction = rep(0,nrow(Assoc_df))

for(i in 1:nrow(Assoc_df)){
  ROW = which(PooledAssoc_df$Trait_A == Assoc_df[i,"Trait_A"] & PooledAssoc_df$Trait_B == Assoc_df[i,"Trait_B"])
  Assoc_df[i,"Direction"] = PooledAssoc_df[ROW,"Direction"][1]
}

Assoc_df$AbsDiff = abs(Assoc_df$Difference) 

p = ggplot(Assoc_df, aes(Trait_A, Trait_B))
p = p + geom_tile(aes(fill = Direction, alpha = AbsDiff), colour = "black")+scale_fill_manual(values = c("#8A0F31", "#36454F","white"))+scale_alpha(range = c(0.3,1))
p = p + xlab("Traits")+ylab("Traits")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
ggsave(paste("TxTAssociation_", Sys.Date(), ".pdf", sep = ""), height = 12, width = 12)

write.csv(PooledAssoc_df, file = paste("AllAsscData_", Sys.Date(), ".csv", sep = ""), row.names = F)
AssociationCounts = data.frame(Trait = character(), Association = character(), Percentage = numeric())
Associations = unique(PooledAssoc_df$Direction)
k = 1
for(i in 1:length(TRAITS)){
  temp = PooledAssoc_df[which(PooledAssoc_df$Trait_A == TRAITS[i] | PooledAssoc_df$Trait_B == TRAITS[i]),]
  for(j in 1:length(Associations)){
    AssociationCounts[k,1] = TRAITS[i]
    AssociationCounts[k,2] = Associations[j]
    AssociationCounts[k,3]= length(which(temp$Direction == Associations[j]))/nrow(temp)
    k = k+1
  }
}
TRAITAORDER = levels(Assoc_df$Trait_A)
AssociationCounts$Trait <- reorder.factor(AssociationCounts$Trait, new.order=TRAITAORDER)
AssociationCounts$Association = factor(AssociationCounts$Association, levels = c("Random","Positive", "Negative"))

a = ggplot(AssociationCounts[order(AssociationCounts$Association, decreasing = T),]
, aes(x = Trait, y = Percentage, fill = Association))+geom_bar(stat = "identity",colour = "black")+scale_fill_manual(values = c("white", "#36454F","#8A0F31"))
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.2,color = "black"))
a = a + xlab("Traits") + ylab("Percent")
ggsave(paste("AsscCounts_", Sys.Date(), ".pdf"))


Random_N=length(which(PooledAssoc_df$Direction == "Random"))/nrow(PooledAssoc_df)
Positive_N=length(which(PooledAssoc_df$Direction == "Positive"))/nrow(PooledAssoc_df)
Negative_N=length(which(PooledAssoc_df$Direction == "Negative"))/nrow(PooledAssoc_df)

Proportion = c(Random_N, Positive_N, Negative_N)
Association = c("Random", "Positive", "Negative")
AssociationCounts = data.frame(Association, Proportion)

a = ggplot(data=AssociationCounts, aes(x=Association, y=Proportion)) + geom_bar(stat="identity")
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
ggsave("AssociationProportions.pdf")

SubsetAssoc_df=PooledAssoc_df[-which(PooledAssoc_df$Direction == "Random"),]
PositiveSub_N=length(which(SubsetAssoc_df$Direction == "Positive"))/nrow(SubsetAssoc_df)
NegativeSub_N=length(which(SubsetAssoc_df$Direction == "Negative"))/nrow(SubsetAssoc_df)

Proportion = c(PositiveSub_N, NegativeSub_N)
Association = c("Positive", "Negative")
SubAssociationCounts = data.frame(Association, Proportion)

a = ggplot(data=SubAssociationCounts, aes(x=Association, y=Proportion)) + geom_bar(stat="identity")
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
ggsave("SubAssociationProportions.pdf")
         
save(PermutedAssoc_list, file, 
			TraitCountPos_m, 
			Assoc_Pos_df,
			TraitCountNeg_m,
			Assoc_Neg_df,
			PooledAssc_df,
			AllAsscData,
			Assoc_m,
			AssociationCounts, 
			SubsetAssoc_df,
			file = "TraitxTrait_Analysis.rda")

