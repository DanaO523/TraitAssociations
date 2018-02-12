############## USED FOR IMPUTING DATA ##############
GrandZero = length(which(TraitData[,1:ncol(TraitData)] == 0))
TotalMatrix = nrow(TraitData)*length(TraitData[,1:ncol(TraitData)])
AvgZero = GrandZero/TotalMatrix

TraitImpute_df = data.frame(Trait = character(), ZeroCount = numeric(), ZeroProp = numeric())
k = 1
for(i in 1:ncol(TraitData)){
  TotalZero = length(which(TraitData[,i] == 0))
  TotalData = length(which(is.na(TraitData[,i]) == FALSE))
  TraitImpute_df[k, 1] = colnames(TraitData)[i]
  TraitImpute_df[k,2] = TotalZero
  TraitImpute_df[k,3] = TotalZero/TotalData
  k = k+1
}

SpeciesImpute_df = data.frame(Species = character(), ZeroCount = numeric(), ZeroProp = numeric())
k = 1
for(i in 1:nrow(TraitData)){
  TotalZero = length(which(TraitData[i,1:ncol(TraitData)] == 0))
  TotalData = length(which(is.na(TraitData[i,1:ncol(TraitData)]) == FALSE))
  SpeciesImpute_df[k, 1] = rownames(TraitData)[i]
  SpeciesImpute_df[k,2] = TotalZero
  SpeciesImpute_df[k,3] = TotalZero/TotalData
  k = k+1
}

SPECIES=rownames(TraitData)
TRAITS = colnames(TraitData)
CxR_ProbMatrix = matrix(0, nrow = length(SPECIES), ncol = length(TRAITS))
colnames(CxR_ProbMatrix) = TRAITS
rownames(CxR_ProbMatrix) = SPECIES

for(i in 1:length(SPECIES)){
  for(j in 1:length(TRAITS)){
    ROW=which(rownames(CxR_ProbMatrix) == SPECIES[i])
    COL = which(colnames(CxR_ProbMatrix) == TRAITS[j])
    SPTab = which(SpeciesImpute_df[,1] == SPECIES[i])
    TRTab = which(TraitImpute_df[,1] == TRAITS[j])
    SPProb = SpeciesImpute_df[SPTab,3]
    TRProb = TraitImpute_df[TRTab,3]
    RxC = SPProb+TRProb
    ProbVal = 1-(AvgZero*RxC)
    CxR_ProbMatrix[ROW,COL] = ProbVal
  }
}

ImputedData_m = matrix(0, nrow = length(SPECIES), ncol = length(TRAITS))
rownames(ImputedData_m) = SPECIES
colnames(ImputedData_m) = TRAITS

for(i in 1:length(SPECIES)){
  for(j in 1:length(TRAITS)){
    IMPROW = which(rownames(ImputedData_m) == SPECIES[i])
    CXRROW = which(rownames(CxR_ProbMatrix) == SPECIES[i])
    IMPCOL = which(colnames(ImputedData_m) == TRAITS[j])
    CXRCOL = which(colnames(CxR_ProbMatrix) == TRAITS[j])
    PROBVALUE = CxR_ProbMatrix[CXRROW, CXRCOL]
    if(PROBVALUE < 0.5){
      ImputedData_m[IMPROW, IMPCOL] = 0
    }else if(PROBVALUE > 0.5){
      ImputedData_m[IMPROW, IMPCOL] = 1
    }else{
      VALUE = sample(c(0,1), size=1, prob=c(PROBVALUE,(1-PROBVALUE)))
      ImputedData_m[IMPROW, IMPCOL] = VALUE
    }
  }
}
write.csv(ImputedData_m, file = "ImputedData.csv")
######## Check Imputation ########

DataCheck = matrix(nrow = nrow(TraitData), ncol = ncol(TraitData))
colnames(DataCheck) = colnames(TraitData)
rownames(DataCheck) = rownames(TraitData)
for(i in 1:nrow(TraitData)){
  for(j in 1:ncol(TraitData)){
    if(is.na(TraitData[i,j]) == F){
      ROW = which(rownames(ImputedData_m) == rownames(TraitData)[i])
      COL = which(colnames(ImputedData_m) == colnames(TraitData)[j])
      DataCheck[i,j] = TraitData[i,j] == ImputedData_m[ROW,COL]
    }else{
      DataCheck[i,j] = NA
    }
  }
}

ColSums = data.frame(Column = numeric(), CountT = numeric(), CountF = numeric())
for(i in 1:ncol(DataCheck)){
  ColSums[i,1] = i
  ColSums[i,2] = length(which(DataCheck[,i] == TRUE))
  ColSums[i,3] = length(which(DataCheck[,i] == FALSE))
}
Total_T=sum(ColSums$CountT)
Total_F=sum(ColSums$CountF)
Prop_T = Total_T/(Total_T+Total_F)

traits = matrix(0, nrow = nrow(TraitData), ncol = ncol(TraitData))
colnames(traits) = colnames(TraitData)
rownames(traits) = rownames(TraitData)

for(i in 1:nrow(traits)){
  for(j in 1:ncol(traits)){
    if(is.na(TraitData[which(rownames(TraitData) == rownames(traits)[i]),which(colnames(TraitData) == colnames(traits)[j])]) == TRUE){
      IMPROW = which(rownames(ImputedData_m) == rownames(traits)[i])
      IMPCOL = which(colnames(ImputedData_m) == colnames(traits)[j])
      traits[i,j] = ImputedData_m[IMPROW,IMPCOL]
    }else{
      traits[i,j] = TraitData[which(rownames(TraitData) == rownames(traits)[i]),which(colnames(TraitData) == colnames(traits)[j])]
    }
  }
}


write.csv(traits, file = paste("TraitData_m_", Sys.Date(), ".csv"))
