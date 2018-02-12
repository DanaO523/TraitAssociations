#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)]  within an isolation environment  #########
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)                                                                                  ###########
#############################################################################################################
# Calculates trait pair counts for permuted data in an isolation environment #
Iso_df = Iso_df[,which(colnames(Iso_df) %in% IsoKeep)]
Iso_df$Species = rownames(Iso_df)
PermAssoc_Iso_TxT_list = list()
for(perm in 1:PERMUTEDVAL){
  PermAssoc_Iso_TxT_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData
  tempDF = merge(tempDF, Iso_df, by = "Species")
  tempDF = tempDF[,-which(colnames(tempDF) == "Species")]
  for(iso in 1:length(ISOLATIONS)){
    PermAssoc_Iso_TxT_list[[perm]][[iso]] = list()
    IsoCOL= which(colnames(tempDF) == ISOLATIONS[[iso]])
    tempComp = tempDF[which(tempDF[,IsoCOL] == 1),which(colnames(tempDF) %in% TRAITS)]
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
        NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
        NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
        tempMatrix_Pos[i,j] = POS
        tempMatrix_Neg[i,j] = NEG1+NEG2
        
      }
    }
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Isolation = ISOLATIONS[[iso]]
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Positive = tempMatrix_Pos
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Negative = tempMatrix_Neg
  }
  print(iso)
}

#save(PermAssoc_Iso_TxT_list, file = "PermAssoc_Iso_TxT_list.rda")

# Calculates trait pair counts for observed data in an isolation environment #
Assoc_Iso_IxT_list = list()
for(iso in 1:length(ISOLATIONS)){
  Assoc_Iso_IxT_list[[iso]] = list()
  IsoCOL= which(colnames(IsoTraits_m) == ISOLATIONS[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_Pos) = colnames(tempComp)
  rownames(tempMatrix_Pos) = colnames(tempComp)
  
  tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_Neg) = colnames(tempComp)
  rownames(tempMatrix_Neg) = colnames(tempComp)
  
  for(i in 1:ncol(tempComp)){
    for(j in 1:ncol(tempComp)){
      POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
      NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
      NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
      tempMatrix_Pos[i,j] = POS
      tempMatrix_Neg[i,j] = NEG1+NEG2
    }
  }
  Assoc_Iso_IxT_list[[iso]]$Isolation = ISOLATIONS[[iso]]
  Assoc_Iso_IxT_list[[iso]]$Positive = tempMatrix_Pos
  Assoc_Iso_IxT_list[[iso]]$Negative = tempMatrix_Neg
  print(iso)
}

###### Isolation Trait pair analysis for Positively Associated traits (+,+) ######
# Determine the number of observed (+/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data #

IsoPermutedList_Pos = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermutedList_Pos[[iso]] = list()
	for(perm in 1:PERMUTEDVAL){
		IsoPermutedList_Pos[[iso]][[perm]] = list()
		IsoPermutedList_Pos[[iso]][[perm]]=PermAssoc_Iso_TxT_list[[perm]][[iso]]$Positive
		print(perm)
	}
	print(iso)
}

TRAITCOLS = melt(IsoPermutedList_Pos[[1]][[1]])

IsoPermMatPos_list = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermMatPos_list[[iso]] = list()
	IsoPermMat = matrix(0, nrow = nrow(TRAITCOLS), ncol = 10003)
	Permuted = c(1:10000)
	COLUMNS = c("Isolation","Trait_A", "Trait_B", Permuted)
	colnames(IsoPermMat) = COLUMNS
	IsoPermMat=data.frame(IsoPermMat)
	traitA = TRAITCOLS[,1]
	traitB = TRAITCOLS[,2]
	isolation = rep(ISOLATIONS[iso], nrow(IsoPermMat))
	IsoPermMat[,1] = isolation
	IsoPermMat[,2] = traitA
	IsoPermMat[,3] = traitB
	for(perm in 1:PERMUTEDVAL){
		test = melt(IsoPermutedList_Pos[[iso]][[perm]])
		IsoPermMat[,(perm+3)] = test[,3]
	}
	IsoPermMatPos_list[[iso]] = IsoPermMat
	print(iso)
}

save(IsoPermMatPos_list, file = "IsoPermMatPos_list.rda")

Iso_TxT_Pos_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
							Observed = numeric(),
							Expected = numeric(),
							Percent = numeric(),
							Count = numeric(),
							Pvalue = numeric())
							
k = 1
for(iso in 1:length(ISOLATIONS)){
  tempData = IsoPermMatPos_list[[iso]]
	for(i in 1:nrow(tempData)){
		isotemp = tempData[i,1]
		trA = tempData[i,2]
		trB = tempData[i,3]
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		tempObs = Assoc_Iso_IxT_list[[iso]]$Positive
		ObsValue = tempObs[which(rownames(tempObs) == trA), which(colnames(tempObs) == trB)] 
		N = sum(temp >= ObsValue)
		p = zapsmall(binconf(N, length(temp), method = "exact"))
		
		Iso_TxT_Pos_df[k,1] = isotemp
		Iso_TxT_Pos_df[k,2] = as.character(trA)
		Iso_TxT_Pos_df[k,3] = as.character(trB)
		Iso_TxT_Pos_df[k,4] = ObsValue
		Iso_TxT_Pos_df[k,5] = mean(temp, na.rm = T)
		Iso_TxT_Pos_df[k,6] = ecdf(temp)(ObsValue)
		Iso_TxT_Pos_df[k,7] = N
		Iso_TxT_Pos_df[k,8] = p[1]
    
		k = k+1
    }
	print(iso)
 }


Iso_TxT_Pos_df$Correction = p.adjust(Iso_TxT_Pos_df$Pvalue, method = "BH")

write.csv(Iso_TxT_Pos_df, file = "Iso_TxT_Pos_df.csv", row.names = F)

IsoPermutedList_Neg = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermutedList_Neg[[iso]] = list()
	for(perm in 1:PERMUTEDVAL){
		IsoPermutedList_Neg[[iso]][[perm]] = list()
		IsoPermutedList_Neg[[iso]][[perm]]=PermAssoc_Iso_TxT_list[[perm]][[iso]]$Negative
		print(perm)
	}
	print(iso)
}

TRAITCOLS = melt(IsoPermutedList_Neg[[1]][[1]])


IsoPermMatNeg_list = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermMatNeg_list[[iso]] = list()
	IsoPermMat = matrix(0, nrow = nrow(TRAITCOLS), ncol = 10003)
	Permuted = c(1:10000)
	COLUMNS = c("Isolation","Trait_A", "Trait_B", Permuted)
	colnames(IsoPermMat) = COLUMNS
	IsoPermMat=data.frame(IsoPermMat)
	traitA = TRAITCOLS[,1]
	traitB = TRAITCOLS[,2]
	isolation = rep(ISOLATIONS[iso], nrow(IsoPermMat))
	IsoPermMat[,1] = isolation
	IsoPermMat[,2] = traitA
	IsoPermMat[,3] = traitB
	for(perm in 1:PERMUTEDVAL){
		test = melt(IsoPermutedList_Neg[[iso]][[perm]])
		IsoPermMat[,(perm+3)] = test[,3]
	}
	IsoPermMatNeg_list[[iso]] = IsoPermMat
	print(iso)
}
save(IsoPermMatNeg_list, file = "IsoPermMatNeg_List.rda")

Iso_TxT_Neg_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
							Observed = numeric(),
							Expected = numeric(),
							Percent = numeric(),
							Count = numeric(),
							Pvalue = numeric())
							
k = 1
for(iso in 1:length(ISOLATIONS)){
  tempData = IsoPermMatNeg_list[[iso]]
	for(i in 1:nrow(tempData)){
		isotemp = tempData[i,1]
		trA = tempData[i,2]
		trB = tempData[i,3]
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		tempObs = Assoc_Iso_IxT_list[[iso]]$Negative
		ObsValue = tempObs[which(rownames(tempObs) == trA), which(colnames(tempObs) == trB)] 
		N = sum(temp >= ObsValue)
		p = zapsmall(binconf(N, length(temp), method = "exact"))
		
		Iso_TxT_Neg_df[k,1] = isotemp
		Iso_TxT_Neg_df[k,2] = as.character(trA)
		Iso_TxT_Neg_df[k,3] = as.character(trB)
		Iso_TxT_Neg_df[k,4] = ObsValue
		Iso_TxT_Neg_df[k,5] = mean(temp, na.rm = T)
		Iso_TxT_Neg_df[k,6] = ecdf(temp)(ObsValue)
		Iso_TxT_Neg_df[k,7] = N
		Iso_TxT_Neg_df[k,8] = p[1]
    
		k = k+1
    }
	print(iso)
 }


Iso_TxT_Neg_df$Correction = p.adjust(Iso_TxT_Neg_df$Pvalue, method = "BH")


write.csv(Iso_TxT_Neg_df, file = "Iso_TxT_Neg_df.csv", row.names = F)

# Summary of permuted data trait pairs within an isolation environment - Postive (+,+) #
Perm_IsoSum_Pos_df = data.frame(Isolation = character(),
                                Trait_A = character(),
                                Trait_B = character(),
                                Mean = numeric(),
                                SD = numeric(),
                                Var = numeric())
k = 1
for(iso in 1:length(ISOLATIONS)){
	tempData = IsoPermMatPos_list[[iso]]
	for(i in 1:nrow(tempData)){
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		Perm_IsoSum_Pos_df[k,1] = as.character(tempData[i,1])
		Perm_IsoSum_Pos_df[k,2] = as.character(tempData[i,2])
		Perm_IsoSum_Pos_df[k,3] = as.character(tempData[i,3])
		Perm_IsoSum_Pos_df[k,4] = mean(temp, na.rm = T)
		Perm_IsoSum_Pos_df[k,5] =sd(temp, na.rm = T)
		Perm_IsoSum_Pos_df[k,6] =var(temp, na.rm = T)
		k = k+1
   }
 }


 ####### Isolation Trait Pair Visualization - Positive (+,+) #######
Iso_TxT_Pos_Summary = merge(Iso_TxT_Pos_df, Perm_IsoSum_Pos_df)
Iso_TxT_Pos_Summary$Difference = Iso_TxT_Pos_Summary$Observed - Iso_TxT_Pos_Summary$Mean
#Iso_TxT_Pos_Summary$Color = rep(0, nrow(Iso_TxT_Pos_Summary))
write.csv(Iso_TxT_Pos_Summary, file = paste("Iso_TxT_Pos_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

Iso_TxT_Stat_Pos_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(ISOLATIONS)){
  tempDF = Iso_TxT_Pos_Summary[which(Iso_TxT_Pos_Summary$Isolation == ISOLATIONS[iso]),]
  Iso_TxT_Stat_Pos_df[iso, 1] = ISOLATIONS[iso]
  Iso_TxT_Stat_Pos_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 4] = mean(tempDF$Observed, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 5] = var(tempDF$Observed, na.rm = TRUE)
}
write.csv(Iso_TxT_Stat_Pos_df, file = paste("Iso_TxT_Stat_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Iso_TxT_Stat_Pos_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("Iso_TxT_Stat_Pos_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Iso_TxT_Stat_Pos_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Iso_TxT_Count_Pos_", Sys.Date(), ".pdf", sep = ""))

###########################################
Perm_IsoSum_Neg_df = data.frame(Isolation = character(),
                                Trait_A = character(),
                                Trait_B = character(),
                                Mean = numeric(),
                                SD = numeric(),
                                Var = numeric())
k = 1
for(iso in 1:length(ISOLATIONS)){
	tempData = IsoPermMatNeg_list[[iso]]
	for(i in 1:nrow(tempData)){
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		Perm_IsoSum_Neg_df[k,1] = as.character(tempData[i,1])
		Perm_IsoSum_Neg_df[k,2] = as.character(tempData[i,2])
		Perm_IsoSum_Neg_df[k,3] = as.character(tempData[i,3])
		Perm_IsoSum_Neg_df[k,4] = mean(temp, na.rm = T)
		Perm_IsoSum_Neg_df[k,5] =sd(temp, na.rm = T)
		Perm_IsoSum_Neg_df[k,6] =var(temp, na.rm = T)
		k = k+1
   }
 }


####### Isolation Trait Pair Visualization - Negitive (+,+) #######
Iso_TxT_Neg_Summary = merge(Iso_TxT_Neg_df, Perm_IsoSum_Neg_df)
Iso_TxT_Neg_Summary$Difference = Iso_TxT_Neg_Summary$Observed - Iso_TxT_Neg_Summary$Mean
#Iso_TxT_Neg_Summary$Color = rep(0, nrow(Iso_TxT_Neg_Summary))
write.csv(Iso_TxT_Neg_Summary, file = paste("Iso_TxT_Neg_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

Iso_TxT_Stat_Neg_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(ISOLATIONS)){
  tempDF = Iso_TxT_Neg_Summary[which(Iso_TxT_Neg_Summary$Isolation == ISOLATIONS[iso]),]
  Iso_TxT_Stat_Neg_df[iso, 1] = ISOLATIONS[iso]
  Iso_TxT_Stat_Neg_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 4] = mean(tempDF$Observed, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 5] = var(tempDF$Observed, na.rm = TRUE)
}
write.csv(Iso_TxT_Stat_Neg_df, file = paste("Iso_TxT_Stat_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Iso_TxT_Stat_Neg_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("Iso_TxT_Stat_Neg_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Iso_TxT_Stat_Neg_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Iso_TxT_Count_Neg_", Sys.Date(), ".pdf", sep = ""))

################################### Detecting Apparent vs Real Association ##################################

### Pick three isolation environments; one with intermediate numbers of species, one with many species and one with few species associated with it.
### Randomly assign a set proportion of traits to positive and negative (separately) and redo 

AvgIso = mean(IsoSums$Counts, na.rm = TRUE)
IsoSp_Sum = IsoSums
IsoSp_Sum$Diff = abs(IsoSp_Sum$Counts - AvgIso)
IntermediateIso = IsoSp_Sum[which(IsoSp_Sum$Diff == min(IsoSp_Sum$Diff)),"Isolation"]
MaxIso = IsoSp_Sum[which(IsoSp_Sum$Counts == max(IsoSp_Sum$Counts)),"Isolation"]
MinIso = IsoSp_Sum[which(IsoSp_Sum$Counts == min(IsoSp_Sum$Counts)),"Isolation"][1]

ISOSUB = c(IntermediateIso, MaxIso, MinIso)
PermAssoc_IsoSub_TxT_list = list()
for(perm in 1:1000){
  PermAssoc_IsoSub_TxT_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData
  tempDF = merge(tempDF, ISODF, by = "Species")
  tempDF = tempDF[,-which(colnames(tempDF) == "Species")]
  for(iso in 1:length(ISOSUB)){
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]] = list()
    IsoCOL= which(colnames(tempDF) == ISOSUB[[iso]])
    tempComp = tempDF[which(tempDF[,IsoCOL] == 1),which(colnames(tempDF) %in% TRAITS)]
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
        NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
        NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
        tempMatrix_Pos[i,j] = POS
        tempMatrix_Neg[i,j] = NEG1+NEG2
        
      }
    }
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Isolation = ISOSUB[[iso]]
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Positive = tempMatrix_Pos
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Negative = tempMatrix_Neg
  }
  print(perm)
}


ChangePercent = c(.75, .5, .25)


Assoc_PosIsoSub_IxT_list = list()
for(iso in 1:length(ISOSUB)){
  Assoc_PosIsoSub_IxT_list[[iso]] = list()
  
  IsoCOL= which(colnames(IsoTraits_m) == ISOSUB[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  for(per in 1:length(ChangePercent)){
    n = ChangePercent[per]*nrow(tempComp)
    ROWCHANGES = sample(nrow(tempComp), n)
    tempDF = tempComp
    Assoc_PosIsoSub_IxT_list[[iso]][[per]] = list()
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        if(i !=j){
          tempDF[ROWCHANGES,c(i,j)] = 1
          POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1))
          tempMatrix_Pos[i,j] = POS  
        }else{
          tempDF[ROWCHANGES,i] = 1
          POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1))
          tempMatrix_Pos[i,j] = POS
        }
        
      }
    }
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Isolation = ISOSUB[[iso]]
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Percentage = ChangePercent[[per]]
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Positive = tempMatrix_Pos
  }
  print(iso)
}
 

Assoc_NegIsoSub_IxT_list = list()
for(iso in 1:length(ISOSUB)){
  Assoc_NegIsoSub_IxT_list[[iso]] = list()
  
  IsoCOL= which(colnames(IsoTraits_m) == ISOSUB[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  for(per in 1:length(ChangePercent)){
    n = ChangePercent[per]*nrow(tempComp)
    ROWCHANGES = sample(nrow(tempComp), n)
    tempDF = tempComp
    Assoc_NegIsoSub_IxT_list[[iso]][[per]] = list()
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        if(i !=j){
          tempDF[ROWCHANGES,i] = 1
          tempDF[ROWCHANGES,j] = 0
          NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
          NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
          tempMatrix_Neg[i,j] = NEG1+NEG2  
        }else{
          tempDF[ROWCHANGES,i] = 1
          NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
          NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
          tempMatrix_Neg[i,j] = NA 
        }
        
      }
    }
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Isolation = ISOSUB[[iso]]
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Percentage = ChangePercent[[per]]
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Negative = tempMatrix_Neg
  }
  print(iso)
}


IsoSub_TxT_Pos_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
                            PerChange = numeric(),
                            Observed = numeric(),
                            Expected = numeric(),
                            Percent = numeric())
k = 1
for(iso in 1:length(ISOSUB)){
  for(per in 1:length(ChangePercent)){
    tempObs = Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Positive
    for(i in 1:length(TRAITS)){
      for(j in 1:length(TRAITS)){
        temp_list = list()
        for(perm in 1:1000){
          VAL = PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Positive[i,j]
          temp_list[perm] = VAL
        }
        templist = unlist(temp_list)
        ObsValue = tempObs[i,j]
        IsoSub_TxT_Pos_df[k,1] = ISOSUB[iso]
        IsoSub_TxT_Pos_df[k,2] = TRAITS[i]
        IsoSub_TxT_Pos_df[k,3] = TRAITS[j]
        IsoSub_TxT_Pos_df[k,4] = Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Percentage
        IsoSub_TxT_Pos_df[k,5] = ObsValue
        IsoSub_TxT_Pos_df[k,6] = mean(templist, na.rm = T)
        IsoSub_TxT_Pos_df[k,7] = ecdf(templist)(ObsValue)
        k = k+1
      }
    }
    print(iso)
  }
}

#############################################################################################################
########### Calculates counts of each type of trait presence/absence [(+/+), (-/-), (+/-), (-/+)]  within an isolation environment  #########
########### for all trait pairs where PP is (+/+), NN is (-/-), PN is (+/-) and NP is (-/+)                                                                                  ###########
#############################################################################################################
# Calculates trait pair counts for permuted data in an isolation environment #
Iso_df = Iso_df[,which(colnames(Iso_df) %in% IsoKeep)]
Iso_df$Species = rownames(Iso_df)
PermAssoc_Iso_TxT_list = list()
for(perm in 1:PERMUTEDVAL){
  PermAssoc_Iso_TxT_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData
  tempDF = merge(tempDF, Iso_df, by = "Species")
  tempDF = tempDF[,-which(colnames(tempDF) == "Species")]
  for(iso in 1:length(ISOLATIONS)){
    PermAssoc_Iso_TxT_list[[perm]][[iso]] = list()
    IsoCOL= which(colnames(tempDF) == ISOLATIONS[[iso]])
    tempComp = tempDF[which(tempDF[,IsoCOL] == 1),which(colnames(tempDF) %in% TRAITS)]
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
        NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
        NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
        tempMatrix_Pos[i,j] = POS
        tempMatrix_Neg[i,j] = NEG1+NEG2
        
      }
    }
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Isolation = ISOLATIONS[[iso]]
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Positive = tempMatrix_Pos
    PermAssoc_Iso_TxT_list[[perm]][[iso]]$Negative = tempMatrix_Neg
  }
  print(iso)
}

#save(PermAssoc_Iso_TxT_list, file = "PermAssoc_Iso_TxT_list.rda")

# Calculates trait pair counts for observed data in an isolation environment #
Assoc_Iso_IxT_list = list()
for(iso in 1:length(ISOLATIONS)){
  Assoc_Iso_IxT_list[[iso]] = list()
  IsoCOL= which(colnames(IsoTraits_m) == ISOLATIONS[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_Pos) = colnames(tempComp)
  rownames(tempMatrix_Pos) = colnames(tempComp)
  
  tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
  colnames(tempMatrix_Neg) = colnames(tempComp)
  rownames(tempMatrix_Neg) = colnames(tempComp)
  
  for(i in 1:ncol(tempComp)){
    for(j in 1:ncol(tempComp)){
      POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
      NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
      NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
      tempMatrix_Pos[i,j] = POS
      tempMatrix_Neg[i,j] = NEG1+NEG2
    }
  }
  Assoc_Iso_IxT_list[[iso]]$Isolation = ISOLATIONS[[iso]]
  Assoc_Iso_IxT_list[[iso]]$Positive = tempMatrix_Pos
  Assoc_Iso_IxT_list[[iso]]$Negative = tempMatrix_Neg
  print(iso)
}

###### Isolation Trait pair analysis for Positively Associated traits (+,+) ######
# Determine the number of observed (+/+)'s for Trait pairs in an isolation environment for the observed data and compare it to the permuted data #

IsoPermutedList_Pos = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermutedList_Pos[[iso]] = list()
	for(perm in 1:PERMUTEDVAL){
		IsoPermutedList_Pos[[iso]][[perm]] = list()
		IsoPermutedList_Pos[[iso]][[perm]]=PermAssoc_Iso_TxT_list[[perm]][[iso]]$Positive
		print(perm)
	}
	print(iso)
}

TRAITCOLS = melt(IsoPermutedList_Pos[[1]][[1]])

IsoPermMatPos_list = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermMatPos_list[[iso]] = list()
	IsoPermMat = matrix(0, nrow = nrow(TRAITCOLS), ncol = 10003)
	Permuted = c(1:10000)
	COLUMNS = c("Isolation","Trait_A", "Trait_B", Permuted)
	colnames(IsoPermMat) = COLUMNS
	IsoPermMat=data.frame(IsoPermMat)
	traitA = TRAITCOLS[,1]
	traitB = TRAITCOLS[,2]
	isolation = rep(ISOLATIONS[iso], nrow(IsoPermMat))
	IsoPermMat[,1] = isolation
	IsoPermMat[,2] = traitA
	IsoPermMat[,3] = traitB
	for(perm in 1:PERMUTEDVAL){
		test = melt(IsoPermutedList_Pos[[iso]][[perm]])
		IsoPermMat[,(perm+3)] = test[,3]
	}
	IsoPermMatPos_list[[iso]] = IsoPermMat
	print(iso)
}

save(IsoPermMatPos_list, file = "IsoPermMatPos_list.rda")

Iso_TxT_Pos_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
							Observed = numeric(),
							Expected = numeric(),
							Percent = numeric(),
							Count = numeric(),
							Pvalue = numeric())
							
k = 1
for(iso in 1:length(ISOLATIONS)){
  tempData = IsoPermMatPos_list[[iso]]
	for(i in 1:nrow(tempData)){
		isotemp = tempData[i,1]
		trA = tempData[i,2]
		trB = tempData[i,3]
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		tempObs = Assoc_Iso_IxT_list[[iso]]$Positive
		ObsValue = tempObs[which(rownames(tempObs) == trA), which(colnames(tempObs) == trB)] 
		N = sum(temp >= ObsValue)
		p = zapsmall(binconf(N, length(temp), method = "exact"))
		
		Iso_TxT_Pos_df[k,1] = isotemp
		Iso_TxT_Pos_df[k,2] = as.character(trA)
		Iso_TxT_Pos_df[k,3] = as.character(trB)
		Iso_TxT_Pos_df[k,4] = ObsValue
		Iso_TxT_Pos_df[k,5] = mean(temp, na.rm = T)
		Iso_TxT_Pos_df[k,6] = ecdf(temp)(ObsValue)
		Iso_TxT_Pos_df[k,7] = N
		Iso_TxT_Pos_df[k,8] = p[1]
    
		k = k+1
    }
	print(iso)
 }


Iso_TxT_Pos_df$Correction = p.adjust(Iso_TxT_Pos_df$Pvalue, method = "BH")

write.csv(Iso_TxT_Pos_df, file = "Iso_TxT_Pos_df.csv", row.names = F)

IsoPermutedList_Neg = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermutedList_Neg[[iso]] = list()
	for(perm in 1:PERMUTEDVAL){
		IsoPermutedList_Neg[[iso]][[perm]] = list()
		IsoPermutedList_Neg[[iso]][[perm]]=PermAssoc_Iso_TxT_list[[perm]][[iso]]$Negative
		print(perm)
	}
	print(iso)
}

TRAITCOLS = melt(IsoPermutedList_Neg[[1]][[1]])


IsoPermMatNeg_list = list()
for(iso in 1:length(ISOLATIONS)){
	IsoPermMatNeg_list[[iso]] = list()
	IsoPermMat = matrix(0, nrow = nrow(TRAITCOLS), ncol = 10003)
	Permuted = c(1:10000)
	COLUMNS = c("Isolation","Trait_A", "Trait_B", Permuted)
	colnames(IsoPermMat) = COLUMNS
	IsoPermMat=data.frame(IsoPermMat)
	traitA = TRAITCOLS[,1]
	traitB = TRAITCOLS[,2]
	isolation = rep(ISOLATIONS[iso], nrow(IsoPermMat))
	IsoPermMat[,1] = isolation
	IsoPermMat[,2] = traitA
	IsoPermMat[,3] = traitB
	for(perm in 1:PERMUTEDVAL){
		test = melt(IsoPermutedList_Neg[[iso]][[perm]])
		IsoPermMat[,(perm+3)] = test[,3]
	}
	IsoPermMatNeg_list[[iso]] = IsoPermMat
	print(iso)
}
save(IsoPermMatNeg_list, file = "IsoPermMatNeg_List.rda")

Iso_TxT_Neg_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
							Observed = numeric(),
							Expected = numeric(),
							Percent = numeric(),
							Count = numeric(),
							Pvalue = numeric())
							
k = 1
for(iso in 1:length(ISOLATIONS)){
  tempData = IsoPermMatNeg_list[[iso]]
	for(i in 1:nrow(tempData)){
		isotemp = tempData[i,1]
		trA = tempData[i,2]
		trB = tempData[i,3]
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		tempObs = Assoc_Iso_IxT_list[[iso]]$Negative
		ObsValue = tempObs[which(rownames(tempObs) == trA), which(colnames(tempObs) == trB)] 
		N = sum(temp >= ObsValue)
		p = zapsmall(binconf(N, length(temp), method = "exact"))
		
		Iso_TxT_Neg_df[k,1] = isotemp
		Iso_TxT_Neg_df[k,2] = as.character(trA)
		Iso_TxT_Neg_df[k,3] = as.character(trB)
		Iso_TxT_Neg_df[k,4] = ObsValue
		Iso_TxT_Neg_df[k,5] = mean(temp, na.rm = T)
		Iso_TxT_Neg_df[k,6] = ecdf(temp)(ObsValue)
		Iso_TxT_Neg_df[k,7] = N
		Iso_TxT_Neg_df[k,8] = p[1]
    
		k = k+1
    }
	print(iso)
 }


Iso_TxT_Neg_df$Correction = p.adjust(Iso_TxT_Neg_df$Pvalue, method = "BH")


write.csv(Iso_TxT_Neg_df, file = "Iso_TxT_Neg_df.csv", row.names = F)

# Summary of permuted data trait pairs within an isolation environment - Postive (+,+) #
Perm_IsoSum_Pos_df = data.frame(Isolation = character(),
                                Trait_A = character(),
                                Trait_B = character(),
                                Mean = numeric(),
                                SD = numeric(),
                                Var = numeric())
k = 1
for(iso in 1:length(ISOLATIONS)){
	tempData = IsoPermMatPos_list[[iso]]
	for(i in 1:nrow(tempData)){
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		Perm_IsoSum_Pos_df[k,1] = as.character(tempData[i,1])
		Perm_IsoSum_Pos_df[k,2] = as.character(tempData[i,2])
		Perm_IsoSum_Pos_df[k,3] = as.character(tempData[i,3])
		Perm_IsoSum_Pos_df[k,4] = mean(temp, na.rm = T)
		Perm_IsoSum_Pos_df[k,5] =sd(temp, na.rm = T)
		Perm_IsoSum_Pos_df[k,6] =var(temp, na.rm = T)
		k = k+1
   }
 }


 ####### Isolation Trait Pair Visualization - Positive (+,+) #######
Iso_TxT_Pos_Summary = merge(Iso_TxT_Pos_df, Perm_IsoSum_Pos_df)
Iso_TxT_Pos_Summary$Difference = Iso_TxT_Pos_Summary$Observed - Iso_TxT_Pos_Summary$Mean
#Iso_TxT_Pos_Summary$Color = rep(0, nrow(Iso_TxT_Pos_Summary))
write.csv(Iso_TxT_Pos_Summary, file = paste("Iso_TxT_Pos_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

Iso_TxT_Stat_Pos_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(ISOLATIONS)){
  tempDF = Iso_TxT_Pos_Summary[which(Iso_TxT_Pos_Summary$Isolation == ISOLATIONS[iso]),]
  Iso_TxT_Stat_Pos_df[iso, 1] = ISOLATIONS[iso]
  Iso_TxT_Stat_Pos_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 4] = mean(tempDF$Observed, na.rm = TRUE)
  Iso_TxT_Stat_Pos_df[iso, 5] = var(tempDF$Observed, na.rm = TRUE)
}
write.csv(Iso_TxT_Stat_Pos_df, file = paste("Iso_TxT_Stat_Pos_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Iso_TxT_Stat_Pos_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("Iso_TxT_Stat_Pos_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Iso_TxT_Stat_Pos_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Iso_TxT_Count_Pos_", Sys.Date(), ".pdf", sep = ""))

###########################################
Perm_IsoSum_Neg_df = data.frame(Isolation = character(),
                                Trait_A = character(),
                                Trait_B = character(),
                                Mean = numeric(),
                                SD = numeric(),
                                Var = numeric())
k = 1
for(iso in 1:length(ISOLATIONS)){
	tempData = IsoPermMatNeg_list[[iso]]
	for(i in 1:nrow(tempData)){
		temp= as.numeric(tempData[i,4:ncol(tempData)])
		Perm_IsoSum_Neg_df[k,1] = as.character(tempData[i,1])
		Perm_IsoSum_Neg_df[k,2] = as.character(tempData[i,2])
		Perm_IsoSum_Neg_df[k,3] = as.character(tempData[i,3])
		Perm_IsoSum_Neg_df[k,4] = mean(temp, na.rm = T)
		Perm_IsoSum_Neg_df[k,5] =sd(temp, na.rm = T)
		Perm_IsoSum_Neg_df[k,6] =var(temp, na.rm = T)
		k = k+1
   }
 }


####### Isolation Trait Pair Visualization - Negitive (+,+) #######
Iso_TxT_Neg_Summary = merge(Iso_TxT_Neg_df, Perm_IsoSum_Neg_df)
Iso_TxT_Neg_Summary$Difference = Iso_TxT_Neg_Summary$Observed - Iso_TxT_Neg_Summary$Mean
#Iso_TxT_Neg_Summary$Color = rep(0, nrow(Iso_TxT_Neg_Summary))
write.csv(Iso_TxT_Neg_Summary, file = paste("Iso_TxT_Neg_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

Iso_TxT_Stat_Neg_df = data.frame(Isolation = character(), AvgDiff = numeric(), VarDiff = numeric(), AvgTraitCounts = numeric(), VarTraitCounts = numeric())
for(iso in 1:length(ISOLATIONS)){
  tempDF = Iso_TxT_Neg_Summary[which(Iso_TxT_Neg_Summary$Isolation == ISOLATIONS[iso]),]
  Iso_TxT_Stat_Neg_df[iso, 1] = ISOLATIONS[iso]
  Iso_TxT_Stat_Neg_df[iso, 2] = mean(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 3] = var(tempDF$Difference, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 4] = mean(tempDF$Observed, na.rm = TRUE)
  Iso_TxT_Stat_Neg_df[iso, 5] = var(tempDF$Observed, na.rm = TRUE)
}
write.csv(Iso_TxT_Stat_Neg_df, file = paste("Iso_TxT_Stat_Neg_df_", Sys.Date(), ".csv", sep = ""), row.names = F)

a = ggplot(Iso_TxT_Stat_Neg_df, aes(x = Isolation, y = AvgDiff))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgDiff-VarDiff, ymax=AvgDiff+VarDiff))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Difference")
ggsave(paste("Iso_TxT_Stat_Neg_", Sys.Date(), ".pdf", sep = ""))

a = ggplot(Iso_TxT_Stat_Neg_df, aes(x = Isolation, y = AvgTraitCounts))+geom_point(pch = 19, fill = "grey60")
a = a + geom_errorbar(aes(ymin=AvgTraitCounts-VarTraitCounts, ymax=AvgTraitCounts+VarTraitCounts))
a = a + geom_hline(yintercept = 0, color = "red", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
a = a + xlab("Isolation") + ylab("Average Trait Counts")
ggsave(paste("Iso_TxT_Count_Neg_", Sys.Date(), ".pdf", sep = ""))

################################### Detecting Apparent vs Real Association ##################################

### Pick three isolation environments; one with intermediate numbers of species, one with many species and one with few species associated with it.
### Randomly assign a set proportion of traits to positive and negative (separately) and redo 

AvgIso = mean(IsoSums$Counts, na.rm = TRUE)
IsoSp_Sum = IsoSums
IsoSp_Sum$Diff = abs(IsoSp_Sum$Counts - AvgIso)
IntermediateIso = IsoSp_Sum[which(IsoSp_Sum$Diff == min(IsoSp_Sum$Diff)),"Isolation"]
MaxIso = IsoSp_Sum[which(IsoSp_Sum$Counts == max(IsoSp_Sum$Counts)),"Isolation"]
MinIso = IsoSp_Sum[which(IsoSp_Sum$Counts == min(IsoSp_Sum$Counts)),"Isolation"][1]

ISOSUB = c(IntermediateIso, MaxIso, MinIso)
PermAssoc_IsoSub_TxT_list = list()
for(perm in 1:1000){
  PermAssoc_IsoSub_TxT_list[[perm]] = list()
  tempDF = PermutedData[[perm]]$CompiledData
  tempDF = merge(tempDF, ISODF, by = "Species")
  tempDF = tempDF[,-which(colnames(tempDF) == "Species")]
  for(iso in 1:length(ISOSUB)){
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]] = list()
    IsoCOL= which(colnames(tempDF) == ISOSUB[[iso]])
    tempComp = tempDF[which(tempDF[,IsoCOL] == 1),which(colnames(tempDF) %in% TRAITS)]
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        POS = length(which(tempComp[,i] == 1 & tempComp[,j] == 1))
        NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
        NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
        tempMatrix_Pos[i,j] = POS
        tempMatrix_Neg[i,j] = NEG1+NEG2
        
      }
    }
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Isolation = ISOSUB[[iso]]
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Positive = tempMatrix_Pos
    PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Negative = tempMatrix_Neg
  }
  print(perm)
}


ChangePercent = c(.75, .5, .25)


Assoc_PosIsoSub_IxT_list = list()
for(iso in 1:length(ISOSUB)){
  Assoc_PosIsoSub_IxT_list[[iso]] = list()
  
  IsoCOL= which(colnames(IsoTraits_m) == ISOSUB[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  for(per in 1:length(ChangePercent)){
    n = ChangePercent[per]*nrow(tempComp)
    ROWCHANGES = sample(nrow(tempComp), n)
    tempDF = tempComp
    Assoc_PosIsoSub_IxT_list[[iso]][[per]] = list()
    
    tempMatrix_Pos = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Pos) = colnames(tempComp)
    rownames(tempMatrix_Pos) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        if(i !=j){
          tempDF[ROWCHANGES,c(i,j)] = 1
          POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1))
          tempMatrix_Pos[i,j] = POS  
        }else{
          tempDF[ROWCHANGES,i] = 1
          POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1))
          tempMatrix_Pos[i,j] = POS
        }
        
      }
    }
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Isolation = ISOSUB[[iso]]
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Percentage = ChangePercent[[per]]
    Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Positive = tempMatrix_Pos
  }
  print(iso)
}
 

Assoc_NegIsoSub_IxT_list = list()
for(iso in 1:length(ISOSUB)){
  Assoc_NegIsoSub_IxT_list[[iso]] = list()
  
  IsoCOL= which(colnames(IsoTraits_m) == ISOSUB[[iso]])
  tempComp = IsoTraits_m[which(IsoTraits_m[,IsoCOL] == 1), which(colnames(IsoTraits_m) %in% TRAITS)]
  
  for(per in 1:length(ChangePercent)){
    n = ChangePercent[per]*nrow(tempComp)
    ROWCHANGES = sample(nrow(tempComp), n)
    tempDF = tempComp
    Assoc_NegIsoSub_IxT_list[[iso]][[per]] = list()
    
    tempMatrix_Neg = matrix(0, ncol = ncol(tempComp), nrow = ncol(tempComp))
    colnames(tempMatrix_Neg) = colnames(tempComp)
    rownames(tempMatrix_Neg) = colnames(tempComp)
    
    for(i in 1:ncol(tempComp)){
      for(j in 1:ncol(tempComp)){
        if(i !=j){
          tempDF[ROWCHANGES,i] = 1
          tempDF[ROWCHANGES,j] = 0
          NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
          NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
          tempMatrix_Neg[i,j] = NEG1+NEG2  
        }else{
          tempDF[ROWCHANGES,i] = 1
          NEG1 = length(which(tempComp[,i] == 1 & tempComp[,j] == 0 ))
          NEG2 = length(which(tempComp[,i] == 0 & tempComp[,j] == 1))        
          tempMatrix_Neg[i,j] = NA 
        }
        
      }
    }
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Isolation = ISOSUB[[iso]]
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Percentage = ChangePercent[[per]]
    Assoc_NegIsoSub_IxT_list[[iso]][[per]]$Negative = tempMatrix_Neg
  }
  print(iso)
}


IsoSub_TxT_Pos_df = data.frame(Isolation = character(),
                            Trait_A = character(),
                            Trait_B = character(),
                            PerChange = numeric(),
                            Observed = numeric(),
                            Expected = numeric(),
                            Percent = numeric())
k = 1
for(iso in 1:length(ISOSUB)){
  for(per in 1:length(ChangePercent)){
    tempObs = Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Positive
    for(i in 1:length(TRAITS)){
      for(j in 1:length(TRAITS)){
        temp_list = list()
        for(perm in 1:1000){
          VAL = PermAssoc_IsoSub_TxT_list[[perm]][[iso]]$Positive[i,j]
          temp_list[perm] = VAL
        }
        templist = unlist(temp_list)
        ObsValue = tempObs[i,j]
        IsoSub_TxT_Pos_df[k,1] = ISOSUB[iso]
        IsoSub_TxT_Pos_df[k,2] = TRAITS[i]
        IsoSub_TxT_Pos_df[k,3] = TRAITS[j]
        IsoSub_TxT_Pos_df[k,4] = Assoc_PosIsoSub_IxT_list[[iso]][[per]]$Percentage
        IsoSub_TxT_Pos_df[k,5] = ObsValue
        IsoSub_TxT_Pos_df[k,6] = mean(templist, na.rm = T)
        IsoSub_TxT_Pos_df[k,7] = ecdf(templist)(ObsValue)
        k = k+1
      }
    }
    print(iso)
  }
}

