ISOLATIONS = colnames(IsolationData)
IsolationCounts_df = data.frame(Isolations = character(), SpeciesCount = numeric())

for(i in 1:length(ISOLATIONS)){
  tempData = IsolationData[,which(colnames(IsolationData) == ISOLATIONS[[i]])]
  IsolationCounts_df[i,1] = ISOLATIONS[i]
  IsolationCounts_df[i,2] = length(which(tempData > 0))
}
#write.csv(IsolationCounts_df, file = paste("IsolationCounts_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

IsoKeep = IsolationCounts_df[which(IsolationCounts_df$SpeciesCount > 3),"Isolations"]
ISOLATIONS = IsoKeep

traits_df = data.frame(traits, check.names = F)
Iso_df = data.frame(IsolationData, check.names = F)

Iso_df = Iso_df[,which(colnames(Iso_df) %in% IsoKeep)]

traits_df$Species = rownames(traits)
Iso_df$Species = rownames(IsolationData)

IsoTraits_m=merge(traits_df, Iso_df, by = "Species")
rownames(IsoTraits_m) = IsoTraits_m[,which(colnames(IsoTraits_m) == "Species")]
IsoTraits_m = IsoTraits_m[,-which(colnames(IsoTraits_m) == "Species")]

IsoDF = Iso_df[,-ncol(Iso_df)]
IsoSums = data.frame(colSums(IsoDF), check.names = F)
IsoSums$Isolation = rownames(IsoSums)
colnames(IsoSums) = c("Counts", "Isolation")


isoP = ggplot(IsoSums, aes(x = Isolation, y = Counts))+geom_bar(stat = "identity", colour = "black", fill = "gray60")
isoP = isoP + xlab("Isolations")+ylab("Species Count")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
ggsave(paste("IsolationCount_", Sys.Date(), ".pdf", sep = ""))


############################################################################
########### Analysis of Isolation Environment/Trait Associations ###########
############################################################################
############## Calculate the number of species isolated from an environment ##############


# Calculate the number of observed [(+/+), (-,-), (+,-), and (-/+)] for each Isolation/Trait pair #
Perm_TxI_Assoc_list = list()
for(perm in 1:PERMUTEDVAL){
  Perm_TxI_Assoc_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData

  tempDF = merge(tempDF, Iso_df, by = "Species")
  tempDF = tempDF[,-which(colnames(tempDF) == "Species")]

  tempMatrix_Pos = matrix(0, ncol = length(TRAITS), nrow = length(ISOLATIONS))
  colnames(tempMatrix_Pos) = TRAITS
  rownames(tempMatrix_Pos) = ISOLATIONS
  
  tempMatrix_Neg = matrix(0, ncol = length(TRAITS), nrow = length(ISOLATIONS))
  colnames(tempMatrix_Neg) = TRAITS
  rownames(tempMatrix_Neg) = ISOLATIONS
  
  for(i in 1:length(ISOLATIONS)){
    for(j in 1:length(TRAITS)){
      ISOCOL = which(colnames(tempDF) == ISOLATIONS[i])
      TRTCOL = which(colnames(tempDF) == TRAITS[j])
      
      POS = length(which(tempDF[,ISOCOL] == 1 & tempDF[,TRTCOL] == 1))
      NEG1 = length(which(tempDF[,ISOCOL] == 1 & tempDF[,TRTCOL] == 0 ))
      NEG2 = length(which(tempDF[,ISOCOL] == 0 & tempDF[,TRTCOL] == 1))
      
      isocol = which(rownames(tempMatrix_Pos) == ISOLATIONS[i])
      trtcol = which(colnames(tempMatrix_Pos) == TRAITS[j])
      
      tempMatrix_Pos[isocol,trtcol] = POS
      tempMatrix_Neg[isocol,trtcol] = NEG1+NEG2
    }
  }
  Perm_TxI_Assoc_list[[perm]]$Positive = tempMatrix_Pos
  Perm_TxI_Assoc_list[[perm]]$Negative = tempMatrix_Neg
  print(perm)
}



###### Determine positive associations between Isolation environments and Traits ######

# Determine the number of observed (+/+)'s for Isolation and Trait for the observed data and compare it to the permuted data #
IsoTrtCountPos_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
rownames(IsoTrtCountPos_m) = ISOLATIONS
colnames(IsoTrtCountPos_m) = TRAITS
for(i in 1:length(ISOLATIONS)){
  for(j in 1:length(TRAITS)){
    ISOCOL = which(colnames(IsoTraits_m) == ISOLATIONS[i])
    TRTCOL = which(colnames(IsoTraits_m) == TRAITS[j])
    POS=length(which(IsoTraits_m[,ISOCOL] == 1 & IsoTraits_m[,TRTCOL] == 1))
    IsoTrtCountPos_m[which(rownames(IsoTrtCountPos_m) == ISOLATIONS[i]),which(colnames(IsoTrtCountPos_m)== TRAITS[j])] = POS
  }
}

ISODF =Iso_df
ISODF$Species = rownames(ISODF)

AssocPos_IxT_df = data.frame(Isolation = character(),
                          Trait = character(),
                          Observed = numeric(),
                          Expected = numeric(),
                          Percent = numeric(),
                          Count = numeric(),
                          Pvalue = numeric())
k = 1
for(i in 1:length(ISOLATIONS)){
  for(j in 1:length(TRAITS)){
    ROWI = which(rownames(IsoTrtCountPos_m) == ISOLATIONS[[i]])
    COLJ= which(colnames(IsoTrtCountPos_m) == TRAITS[[j]])
    
    tempObs = IsoTrtCountPos_m[i,j]
    tempList = list()
    for(perm in 1:PERMUTEDVAL){
      tempDF1 = Perm_TxI_Assoc_list[[perm]]$Positive
      OBS = tempDF1[ROWI,COLJ]
      tempList[perm] = OBS
    }
    temp = unlist(tempList)
    N = sum(temp >= tempObs)
    p = zapsmall(binconf(N, length(temp), method = "exact"))
    
    AssocPos_IxT_df[k, 1] = ISOLATIONS[[i]]
    AssocPos_IxT_df[k, 2] = TRAITS[[j]]
    AssocPos_IxT_df[k, 3] = tempObs
    AssocPos_IxT_df[k, 4] = mean(temp, na.rm = T)
    AssocPos_IxT_df[k, 5] = ecdf(temp)(tempObs)
    AssocPos_IxT_df[k, 6] = N
    AssocPos_IxT_df[k, 7] = p[1]
    k = k + 1
  }
}

AssocPos_IxT_df$Sig = AssocPos_IxT_df$Pvalue <= 0.05
AssocPos_IxT_df$Correction = p.adjust(AssocPos_IxT_df$Pvalue, method = "BH")

# Determine which traits are positive (+,+) significant in an environment #

SigTraits_IxT_Pos_list = list()
for(i in 1:length(ISOLATIONS)){
  SigTraits_IxT_Pos_list[[i]] = list()
  temp_df = AssocPos_IxT_df[which(AssocPos_IxT_df$Isolation == ISOLATIONS[i]),]
  SigTraits_IxT_Pos_list[[i]]$Isolation = ISOLATIONS[i]
  SigTraits_IxT_Pos_list[[i]]$Traits = temp_df[which(temp_df$Sig == "TRUE"),"Trait"]
}

AssocPos_IxT_df$TraitDiff = AssocPos_IxT_df[,3] - AssocPos_IxT_df[,4]
#Assoc_IxT_Pos_df$ID = c(1:nrow(Assoc_IxT_Pos_df))
Assoc_IxT_Pos_mean = mean(AssocPos_IxT_df$TraitDiff, na.rm = T)

###### Determine negative associations between Isolation environments and Traits ######
# Determine the number of observed (+/-)'s for Isolation and Trait for the observed data and compare it to the permuted data #
IsoTrtCountNeg_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
rownames(IsoTrtCountNeg_m) = ISOLATIONS
colnames(IsoTrtCountNeg_m) = TRAITS
for(i in 1:length(ISOLATIONS)){
  for(j in 1:length(TRAITS)){
    ISOCOL = which(colnames(IsoTraits_m) == ISOLATIONS[i])
    TRTCOL = which(colnames(IsoTraits_m) == TRAITS[j])
    NEG1=length(which(IsoTraits_m[,ISOCOL] == 0 & IsoTraits_m[,TRTCOL] == 1)) 
    NEG2 = length(which(IsoTraits_m[,ISOCOL] == 1 & IsoTraits_m[,TRTCOL] == 0))
    NEG = NEG1+NEG2
    IsoTrtCountNeg_m[which(rownames(IsoTrtCountNeg_m) == ISOLATIONS[i]),which(colnames(IsoTrtCountNeg_m)== TRAITS[j])] = NEG
  }
}

AssocNeg_IxT_df = data.frame(Isolation = character(),
                             Trait = character(),
                             Observed = numeric(),
                             Expected = numeric(),
                             Percent = numeric(),
                             Count = numeric(),
                             Pvalue = numeric())
k = 1
for(i in 1:length(ISOLATIONS)){
  for(j in 1:length(TRAITS)){
    ROWI = which(rownames(IsoTrtCountNeg_m) == ISOLATIONS[[i]])
    COLJ= which(colnames(IsoTrtCountNeg_m) == TRAITS[[j]])
    
    tempObs = IsoTrtCountNeg_m[i,j]
    tempList = list()
    for(perm in 1:PERMUTEDVAL){
      tempDF1 = Perm_TxI_Assoc_list[[perm]]$Negative
      ROWNUMB = which(rownames(tempDF1) == ISOLATIONS[[i]])
      COLNUMB = which(colnames(tempDF1) == TRAITS[[j]])
      OBS = tempDF1[ROWI,COLJ]
      tempList[perm] = OBS
    }
    temp = unlist(tempList)
    N = sum(temp >= tempObs)
    p = zapsmall(binconf(x = N, n = length(temp), alpha = 0.05,method = "exact"))
    
    AssocNeg_IxT_df[k, 1] = ISOLATIONS[[i]]
    AssocNeg_IxT_df[k, 2] = TRAITS[[j]]
    AssocNeg_IxT_df[k, 3] = tempObs
    AssocNeg_IxT_df[k, 4] = mean(temp, na.rm = T)
    AssocNeg_IxT_df[k, 5] = ecdf(temp)(tempObs)
    AssocNeg_IxT_df[k, 6] = N
    AssocNeg_IxT_df[k, 7] = p[1]
    k = k + 1
  }
}

AssocNeg_IxT_df$Sig = AssocNeg_IxT_df$Pvalue <= 0.05
AssocPos_IxT_df$Correction = p.adjust(AssocPos_IxT_df$Pvalue, method = "bonferroni")

# Determine which traits are positive (+,+) significant in an environment #
SigTraits_IxT_Neg_list = list()
for(i in 1:length(ISOLATIONS)){
  SigTraits_IxT_Neg_list[[i]] = list()
  temp_df = AssocNeg_IxT_df[which(AssocNeg_IxT_df$Isolation == ISOLATIONS[i]),]
  SigTraits_IxT_Neg_list[[i]]$Isolation = ISOLATIONS[i]
  SigTraits_IxT_Neg_list[[i]]$Traits = temp_df[which(temp_df$Sig == "TRUE"),"Trait"]
}

AssocNeg_IxT_df$TraitDiff = AssocNeg_IxT_df[,3] - AssocNeg_IxT_df[,4]
#Assoc_IxT_Neg_df$ID = c(1:nrow(Assoc_IxT_Neg_df))
Assoc_IxT_Neg_mean = mean(AssocNeg_IxT_df$TraitDiff, na.rm = T)


AssocPos_IxT_df$Direction = rep("Random", nrow(AssocPos_IxT_df))
AssocNeg_IxT_df$Direction = rep("Random", nrow(AssocNeg_IxT_df))
Assoc_IxT_Pos_df = AssocPos_IxT_df
Assoc_IxT_Neg_df = AssocNeg_IxT_df

for(i in 1:nrow(Assoc_IxT_Pos_df)){
  if(Assoc_IxT_Pos_df[i,"Sig"] == TRUE)
    Assoc_IxT_Pos_df[i,"Direction"] = "Positive"
}

for(i in 1:nrow(Assoc_IxT_Neg_df)){
  if(Assoc_IxT_Neg_df[i,"Sig"] == TRUE)
    Assoc_IxT_Neg_df[i,"Direction"] = "Negative"
}

Assoc_IxT_df = data.frame(rbind(Assoc_IxT_Pos_df, Assoc_IxT_Neg_df))
Assoc_IxT_Data = matrix(0, nrow = 1, ncol = ncol(Assoc_IxT_df))
colnames(Assoc_IxT_Data) = colnames(Assoc_IxT_df)
Assoc_IxT_Data = data.frame(Assoc_IxT_Data)

k = 1

for(i in 1:length(ISOLATIONS)){
  for(j in 1:length(TRAITS)){
    tempData = Assoc_IxT_df[which(Assoc_IxT_df$Isolation == ISOLATIONS[i] & Assoc_IxT_df$Trait == TRAITS[j]),]
    if(length(unique(tempData$Direction)) == 1){
      Assoc_IxT_Data[k,] = tempData[1,]
      k = k + 1
    }else if(length(unique(tempData$Direction)) >1 & length(which(tempData$Direction == "Random")) == 1){
      Assoc_IxT_Data[k,] = tempData[which(tempData$Direction != "Random"),]
      k = k + 1
    }else if(length(unique(tempData$Direction)) >1 & length(which(tempData$Direction == "Random")) > 1){
      Assoc_IxT_Data[k, ] = tempData[1,]
      Assoc_IxT_Data[k,"Direction"] = "Shit"
      k = k + 1
      Assoc_IxT_Data[k, ] = tempData[2,]
      Assoc_IxT_Data[k,"Direction"] = "Shit"
      k = k+1
    }
  }
}


Assoc_IxT_Data$AbsDiff = abs(Assoc_IxT_Data$TraitDiff)
write.csv(Assoc_IxT_Data, file = "Assoc_IxT_df.csv", row.names = F)

NEGROW = which(Assc_IxT_Data$Direction == "Negative")
RANDROW = which(Assc_IxT_Data$Direction == "Random")

test  = Assc_IxT_Data
for(i in 1:length(NEGROW)){
  test[NEGROW[[i]],"AbsDiff"] = (test[NEGROW[[i]],"AbsDiff"])*-1
}

for(i in 1:length(RANDROW)){
  test[RANDROW[[i]],"AbsDiff"] = (test[RANDROW[[i]],"AbsDiff"])*0
}

IxT_Assoc_m = matrix(0, nrow = length(TRAITS), ncol = length(ISOLATIONS))
colnames(IxT_Assoc_m) = ISOLATIONS
rownames(IxT_Assoc_m) = TRAITS

for(i in 1:nrow(test)){
  IxT_Assoc_m[which(rownames(IxT_Assoc_m) == test[i,"Trait"]),which(colnames(IxT_Assoc_m) == test[i,"Isolation"])] = test[i,"AbsDiff"]
}

IxT_fit = pvclust(IxT_Assoc_m, method.dist = "euclidean", method.hclust = "ward.D", nboot = 1000, store = T)
IxT_clusters= pvpick(IxT_fit, max.only = F)
pdf("IsolationClusters.pdf")
plot(IxT_fit, print.num = F)
pvrect(IxT_fit, max.only = F)
dev.off()

Assoc_IxT_Counts = data.frame(Isolation = character(), Association = character(), Percentage = numeric())
Associations = unique(Assc_IxT_Data$Direction)

k = 1
for(i in 1:length(ISOLATIONS)){
  temp = Assc_IxT_Data[which(Assc_IxT_Data$Isolation == ISOLATIONS[i]),]
  for(j in 1:length(Associations)){
    Assoc_IxT_Counts[k,1] = ISOLATIONS[i]
    Assoc_IxT_Counts[k,2] = Associations[j]
    Assoc_IxT_Counts[k,3]= length(which(temp$Direction == Associations[j]))/nrow(temp)
    k = k+1
  }
}
ISOORDER = levels(IxT_Assoc_df$Isolation)
Assoc_IxT_Counts$Isolation <- reorder.factor(Assoc_IxT_Counts$Isolation, new.order=ISOORDER)
Assoc_IxT_Counts$Association = factor(Assoc_IxT_Counts$Association, levels = c("Random","Positive", "Negative"))

a = ggplot(Assoc_IxT_Counts[order(Assoc_IxT_Counts$Association, decreasing = T),]
           , aes(x = Isolation, y = Percentage, fill = Association))+geom_bar(stat = "identity",colour = "black")+scale_fill_manual(values = c("white", "#36454F","#8A0F31"))
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.2,color = "black"))
a = a + xlab("Isolation") + ylab("Percent")
ggsave(paste("IxT_AsscCounts_", Sys.Date(), ".pdf"))


TxI_fit = pvclust(t(IxT_Assoc_m), method.dist = "euclidean", method.hclust = "ward.D", nboot = 1000, store = F)
plot(TxI_fit, print.num = F)
pvrect(TxI_fit, max.only = T)
TxI_clusters= pvpick(TxI_fit, max.only = F)

ISODROP = Assoc_IxT_Counts[which(Assoc_IxT_Counts$Association == "Random" & Assoc_IxT_Counts$Percentage == 1),"Isolation"]
IxT_AssocSub_m = IxT_Assoc_m[,-which(colnames(IxT_Assoc_m) %in% ISODROP)]

Tr_hclust =hclust( dist(IxT_AssocSub_m, method = "euclidean"), method = "ward.D" )
Iso_hclust =hclust( dist(t(IxT_AssocSub_m), method = "euclidean"), method = "ward.D")

T_ord <- Tr_hclust$order
I_ord <- Iso_hclust$order
IxT_pd <- as.data.frame(IxT_AssocSub_m)
IxT_pd$Trait <- rownames(IxT_AssocSub_m)
IxT_pd.m <- melt( IxT_pd, id.vars = "Trait" )
colnames(IxT_pd.m) = c("Trait", "Isolation", "Difference")

IxT_pd.m$Trait <- factor(IxT_pd.m$Trait, levels = rownames(IxT_AssocSub_m)[T_ord])
IxT_pd.m$Isolation <- factor(IxT_pd.m$Isolation, levels = colnames(IxT_AssocSub_m)[I_ord])

# Melt the correlation matrix
IxT_Assoc_df = IxT_pd.m
IxT_Assoc_df$Direction = rep(0,nrow(IxT_Assoc_df))

for(i in 1:nrow(IxT_Assoc_df)){
  ROW = which(Assc_IxT_Data$Trait == IxT_Assoc_df[i,"Trait"] & Assc_IxT_Data$Isolation == IxT_Assoc_df[i,"Isolation"])
  IxT_Assoc_df[i,"Direction"] = Assc_IxT_Data[ROW,"Direction"][1]
}

IxT_Assoc_df$AbsDiff = abs(IxT_Assoc_df$Difference) 

p = ggplot(IxT_Assoc_df, aes(Isolation, Trait))
p = p + geom_tile(aes(fill = Direction, alpha = AbsDiff), colour = "black")+scale_fill_manual(values = c("#8A0F31", "#36454F","white"))+scale_alpha(range = c(0.4,1))
p = p + xlab("Isolation")+ylab("Traits")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
ggsave(paste("IxTAssociation_", Sys.Date(), ".pdf", sep = ""), height = 12, width = 12)
