options(stringsAsFactors = FALSE)

library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("Hmisc")
library("reshape")

IsoTraits_m = read.csv("IsoTrait_m_2017-01-30.csv", header = TRUE, check.names = F)
load("IsoPermMatPos_list.rda")
load("TxI_AssocAnalysis_2017-01-26.rda")
load("PermutedData.rda")

AvgIso = mean(IsoSums$Counts, na.rm = TRUE)
IsoSp_Sum = IsoSums
IsoSp_Sum$Diff = abs(IsoSp_Sum$Counts - AvgIso)
IntermediateIso = IsoSp_Sum[which(IsoSp_Sum$Diff == min(IsoSp_Sum$Diff)),"Isolation"]
MaxIso = IsoSp_Sum[which(IsoSp_Sum$Counts == max(IsoSp_Sum$Counts)),"Isolation"]
MinIso = IsoSp_Sum[which(IsoSp_Sum$Counts == min(IsoSp_Sum$Counts)),"Isolation"][1]

ISOSUB = c(IntermediateIso, MaxIso, MinIso, "Clinical", "Beer", "Plant", "Beetle")

ISO2USE = c(14,7,19, 2, 44, 36, 18)

PermMatKeep = list()

for(iso in 1:length(ISO2USE)){
	PermMatKeep[[iso]] = IsoPermMatPos_list[[ISO2USE[[iso]]]]
}
			
RandomAssoc_PosIsoSub_IxT_list = list()
for(iso in 1:length(ISOSUB)){
	RandomAssoc_PosIsoSub_IxT_list[[iso]] = list()
	for(perm in 1:1000){
		RandomAssoc_PosIsoSub_IxT_list[[iso]][[perm]] = list()

		COLS = which(colnames(IsoTraits_m) == ISOSUB[[iso]])
		n =length(which(IsoTraits_m[ ,COLS] == 1))
#		ROWKEEP = which(IsoTraits_m[ ,COLS] == 1)
		ROWCHANGES = sample(nrow(IsoTraits_m), n)

		tempDF = IsoTraits_m[ROWCHANGES,1:48]
	
 
		tempMatrix_Pos = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
		colnames(tempMatrix_Pos) = colnames(tempDF)
		rownames(tempMatrix_Pos) = colnames(tempDF)
    
		for(i in 1:ncol(tempDF)){
			for(j in 1:ncol(tempDF)){
				POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1 ))
				tempMatrix_Pos[i,j] = POS
			}
		}
		RandomAssoc_PosIsoSub_IxT_list[[iso]][[perm]]$Isolation = ISOSUB[[iso]]
		RandomAssoc_PosIsoSub_IxT_list[[iso]][[perm]]$Positive = tempMatrix_Pos
		print(iso)
	}
}

TRAITCOLS2 = PermMatKeep[[1]][,2:3]

RandomIsoTrait_Pos_list = list()
for(iso in 1:length(ISOSUB)){
	RandomIsoTrait_Pos_list[[iso]] = list()
	random_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 1002)
	Permuted = c(1:1000)
	COLUMNS = c("Trait_A", "Trait_B", Permuted)
	colnames(random_m) = COLUMNS
	random_m=data.frame(random_m)
	traitA = TRAITCOLS2[,1]
	traitB = TRAITCOLS2[,2]
	random_m[,1] = traitA
	random_m[,2] = traitB
	for(perm in 1:1000){
		test = RandomAssoc_PosIsoSub_IxT_list[[iso]][[perm]]$Positive
		test = melt(test)
		random_m[,(perm+2)] = test[,3]
	}
		RandomIsoTrait_Pos_list[[iso]]$Isolation = ISOSUB[[iso]]
		RandomIsoTrait_Pos_list[[iso]]$Matrix = random_m
	print(iso)
}


RandomExpected_Pos_list = list()
for(iso in 1:length(ISOSUB)){ 
	RandomExpected_Pos_list[[iso]] = list()
	tempData = PermMatKeep[[iso]]
	random_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 3)
	COLUMNS = c("Trait_A", "Trait_B", "Expected")
	colnames(random_m) = COLUMNS
	random_m=data.frame(random_m)
	traitA = tempData[,2]
	traitB = tempData[,3]
	random_m[,1] = traitA
	random_m[,2] = traitB
	random_m[,3] = apply(tempData[4:ncol(tempData)],1,mean)
	RandomExpected_Pos_list[[iso]]$Isolation = ISOSUB[[iso]]
	RandomExpected_Pos_list[[iso]]$Data.Frame = random_m
	print(iso)
}

RandomDiff_Pos_list = list()
for(iso in 1:length(ISOSUB)){
	RandomDiff_Pos_list[[iso]] = list()
	tempData = RandomIsoTrait_Pos_list[[iso]]$Matrix
	test =RandomExpected_Pos_list[[iso]]$Data.Frame
	random_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 2)
	COLUMNS = c("Trait_A", "Trait_B")
	colnames(random_m) = COLUMNS
	random_m=data.frame(random_m)

	traitA = tempData[,1]
	traitB = tempData[,2]
	random_m[,1] = traitA
	random_m[,2] = traitB
	WOOT = tempData[,3:ncol(tempData)]-test[,3]
	random_m = cbind(random_m, WOOT)
	RandomDiff_Pos_list[[iso]]$Isolation = ISOSUB[[iso]]
	RandomDiff_Pos_list[[iso]]$Matrix = random_m

}

AvgDiffPos_df = matrix(0, ncol = 1000, nrow = length(ISOSUB))
for(iso in 1:length(ISOSUB)){
	tempDF = RandomDiff_Pos_list[[iso]]$Matrix
	YEAH = t(data.frame(apply(tempDF[,3:ncol(tempDF)], 2, mean)))
	AvgDiffPos_df[iso,] = YEAH
}

AvgDiffPos_df = data.frame(AvgDiffPos_df)
AvgDiffPos_df = cbind(data.frame(ISOSUB), AvgDiffPos_df)

AvgDiffStat_Pos_df = data.frame(Isolation=AvgDiffPos_df[,1], 
													AvgDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, mean),
													VarDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, var),
													StErrDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, function(x) sd(x, na.rm=TRUE) /sqrt(length(x[!is.na(x)]))))
													

write.csv(AvgDiffStat_Pos_df, file = paste("Random_Iso_TxT_Pos_Check_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

save(PermMatKeep, 
			RandomAssoc_PosIsoSub_IxT_list,
			RandomIsoTrait_Pos_list,
			RandomExpected_Pos_list,
			RandomDiff_Pos_list,
			AvgDiffPos_df,
			AvgDiffStat_Pos_df,
			file = "Random_Pos_Iso_TxT_Check.rda")

			
NumbtoSample = c(min(IsoSp_Sum$Counts):max(IsoSp_Sum$Counts))

SamplePull_PosI_list = list()
for(samp in 1:length(NumbtoSample)){
	SamplePull_PosI_list[[samp]] = list()
	for(perm in 1:1000){
		SamplePull_PosI_list[[samp]][[perm]] = list()

		n =NumbtoSample[[samp]]
		ROWCHANGES = sample(nrow(IsoTraits_m), n)

		tempDF = IsoTraits_m[ROWCHANGES,1:48]
	
 
		tempMatrix_Pos = matrix(0, ncol = ncol(tempDF), nrow = ncol(tempDF))
		colnames(tempMatrix_Pos) = colnames(tempDF)
		rownames(tempMatrix_Pos) = colnames(tempDF)
    
		for(i in 1:ncol(tempDF)){
			for(j in 1:ncol(tempDF)){
				POS = length(which(tempDF[,i] == 1 & tempDF[,j] == 1 ))
				tempMatrix_Pos[i,j] = POS
			}
		}
		SamplePull_PosI_list[[samp]][[perm]]$Isolation = NumbtoSample[[samp]]
		SamplePull_PosI_list[[samp]][[perm]]$Positive = tempMatrix_Pos
	}
	print(samp)
}

TRAITCOLS2 = PermMatKeep[[1]][,2:3]
RandomSamplePull_Pos_list = list()
for(samp in 1:length(NumbtoSample)){
	RandomSamplePull_Pos_list[[samp]] = list()
	random_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 1002)
	Permuted = c(1:1000)
	COLUMNS = c("Trait_A", "Trait_B", Permuted)
	colnames(random_m) = COLUMNS
	random_m=data.frame(random_m)
	traitA = TRAITCOLS2[,1]
	traitB = TRAITCOLS2[,2]
	random_m[,1] = traitA
	random_m[,2] = traitB
	for(perm in 1:1000){
		test = melt(SamplePull_PosI_list[[samp]][[perm]]$Positive)
		random_m[,(perm+2)] = test[,3]
	}
		RandomSamplePull_Pos_list[[samp]]$Sampling = NumbtoSample[[samp]]
		RandomSamplePull_Pos_list[[samp]]$Matrix = random_m
	print(samp)
}

PERMUTEDVAL = 10000
TRAITS = PermMatKeep[[iso]][,2]
TRAITS = unique(TRAITS)

PermutedAssoc_list = list()
for(numb in 1:length(NumbtoSample)){
	PermutedAssoc_list[[numb]] = list()
	for(perm in 1:PERMUTEDVAL){
		PermutedAssoc_list[[numb]][[perm]] = list()
		tempDF = PermutedData[[perm]]$CompiledData
		ROWS= sample(1:nrow(tempDF), size = NumbtoSample[[numb]])
		tempDF = tempDF[,-which(colnames(tempDF) == "Species")]
		tempDF=tempDF[ROWS,which(colnames(tempDF) %in% TRAITS)]
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
		PermutedAssoc_list[[numb]][[perm]]$Number = NumbtoSample[[numb]]
		PermutedAssoc_list[[numb]][[perm]]$Positive = tempMatrix_Pos
		PermutedAssoc_list[[numb]][[perm]]$Negative = tempMatrix_Neg
	}
	print(NumbtoSample[[numb]])
}

save(PermutedAssoc_list, file = paste("SampledSubset_PermutedAssoc_list_", Sys.Date(), ".rda", sep = ""))

TRAITCOLS = melt(PermutedAssoc_list[[1]]$Positive)

RandomPermutedMatrix_list = list()
for(numb in 1:length(NumbtoSample)){
	RandomPermutedMatrix_list[[numb]] = list()
	for(perm in 1:1000)
		RandomPermutedMatrix_list[[numb]][[perm]] = list()
		PermPos_m = matrix(0, nrow = nrow(TRAITCOLS), ncol = 10002)
		Permuted = c(1:10000)
		COLUMNS = c("Trait_A", "Trait_B", Permuted)
		colnames(PermPos_m) = COLUMNS
		PermPos_m=data.frame(PermPos_m)
		PermPos_m[,1] = TRAITCOLS[,1]
		PermPos_m[,2] = TRAITCOLS[,2]
			for(perm in 1:PERMUTEDVAL){
				test = melt(PermutedAssoc_list[[perm]]$Positive)
				PermPos_m[,(perm+2)] = test[,3]
			}


PermPosExp_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 3)
COLUMNS = c("Trait_A", "Trait_B", "Expected")
colnames(PermPosExp_m) = COLUMNS
PermPosExp_m=data.frame(PermPosExp_m)
PermPosExp_m[,1] = TRAITCOLS[,1]
PermPosExp_m[,2] = TRAITCOLS[,2]
PermPosExp_m[,3] = apply(tempData[4:ncol(tempData)],1,mean)


RandomSampleDiff_Pos_list = list()
for(samp in 1:length(NumbtoSample)){
	RandomDiff_Pos_list[[samp]] = list()
	tempData = RandomSamplePull_Pos_list[[samp]]$Matrix
	test =PermPosExp_m[[samp]]$Data.Frame
	random_m = matrix(0, nrow = nrow(TRAITCOLS2), ncol = 2)
	COLUMNS = c("Trait_A", "Trait_B")
	colnames(random_m) = COLUMNS
	random_m=data.frame(random_m)

	traitA = tempData[,1]
	traitB = tempData[,2]
	random_m[,1] = traitA
	random_m[,2] = traitB
	WOOT = tempData[,3:ncol(tempData)]-test[,3]
	random_m = cbind(random_m, WOOT)
	RandomDiff_Pos_list[[samp]]$Isolation = NumbtoSample[[samp]]
	RandomDiff_Pos_list[[samp]]$Matrix = random_m

}

AvgDiffPos_df = matrix(0, ncol = 1000, nrow = length(NumbtoSample))
for(samp in 1:length(NumbtoSample)){
	tempDF = RandomDiff_Pos_list[[samp]]$Matrix
	YEAH = t(data.frame(apply(tempDF[,3:ncol(tempDF)], 2, mean)))
	AvgDiffPos_df[samp,] = YEAH
}

AvgDiffPos_df = data.frame(AvgDiffPos_df)
AvgDiffPos_df = cbind(data.frame(NumbtoSample), AvgDiffPos_df)

AvgDiffStat_Pos_df = data.frame(Isolation=AvgDiffPos_df[,1], 
													AvgDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, mean),
													VarDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, var),
													StErrDiff = apply(AvgDiffPos_df[,2:ncol(AvgDiffPos_df)], 1, function(x) sd(x, na.rm=TRUE) /sqrt(length(x[!is.na(x)]))))
													

write.csv(AvgDiffStat_Pos_df, file = paste("Random_Iso_TxT_Pos_Check_Summary_", Sys.Date(), ".csv", sep = ""), row.names = F)

save(PermMatKeep, 
			RandomAssoc_PosIsoSub_IxT_list,
			RandomIsoTrait_Pos_list,
			RandomExpected_Pos_list,
			RandomDiff_Pos_list,
			AvgDiffPos_df,
			AvgDiffStat_Pos_df,
			file = "Random_Pos_Iso_TxT_Check.rda")



