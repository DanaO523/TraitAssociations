options(stringsAsFactors = FALSE)
library("grofit")
library("doBy")
library("RColorBrewer")
library("pvclust")
library("ggplot2")
library("fBasics")
library("ape")
library("pheatmap")
library("reshape")

Imputed = read.csv("AllAssocImpute.csv", header = TRUE)
all1 = read.csv("AllAssocAll1.csv", header = TRUE)
all0 = read.csv("AllAssocAll0.csv", header = TRUE)
RawData = read.csv("RawData.csv", header = TRUE)
AllImputed = read.csv("ImputedData.csv", header = TRUE)

GrowthQual_m1 = read.csv("GrowthQual_m1.csv", header = TRUE, check.names = F)
GrowthQual_m0 = read.csv("GrowthQual_m0.csv", header = TRUE, check.names = F)

qual1_1 = length(which(GrowthQual_m1 == 1 | GrowthQual_m1 == 0.5))
qual1_0 = length(which(GrowthQual_m1 == 0))

qual0_1 = length(which(GrowthQual_m0 == 1| GrowthQual_m1 == 0.5))
qual0_0 = length(which(GrowthQual_m0 == 0))

traits_1 = length(which(AllImputed == 1))
traits_0 = length(which(AllImputed == 0))

ImputedData_df = data.frame(matrix(0,ncol = 3, nrow = 6))
colnames(ImputedData_df) = c("Data", "Values", "Counts")
ImputedData_df$Data = c("Qual1", "Qual1", "Qual0", "Qual0", "Imputed", "Imputed")
ImputedData_df$Values = c(1, 0 , 1, 0, 1, 0)
ImputedData_df$Counts = c(qual1_1, qual1_0 , qual0_1, qual0_0, traits_1, traits_0)
totalCells = nrow(AllImputed)*(ncol(AllImputed)-1)/2


a = ggplot(data=ImputedData_df, aes(x=Data, y=Counts, fill=factor(Values))) +geom_bar(stat="identity")+scale_fill_manual(values = c("#8c510a", "#01665e"))
a = a + geom_hline(yintercept = totalCells, color = "#dfc27d", linetype  = 2)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
ggsave("ImputationComparisons.pdf")

ImpComp=data.frame(Trait_A = character(), Trait_B = character(), AllEqual = character(), Agree0 = character(), Agree1 = character(),ImputDir = character()) 
k = 1
for(i in 1:nrow(Imputed)){
  temp1Dir=all1[which(all1$Trait_A == Imputed[i,"Trait_A"] & all1$Trait_B == Imputed[i,"Trait_B"]),"Direction"]
  temp0Dir=all0[which(all0$Trait_A == Imputed[i,"Trait_A"] & all0$Trait_B == Imputed[i,"Trait_B"]),"Direction"]
  tempImpDir = Imputed[i,"Direction"]
  if(length(temp1Dir) >= length(temp0Dir)){
    for(j in 1:length(temp1Dir)){
      for(l in 1:length(temp0Dir)){
        ImpComp[k,1] = Imputed[i,"Trait_A"]
        ImpComp[k,2] = Imputed[i,"Trait_B"]
        ImpComp[k,3] = (tempImpDir == temp0Dir[l]) == (tempImpDir == temp1Dir[j])
        ImpComp[k,4] = (tempImpDir == temp0Dir[l])
        ImpComp[k,5] = (tempImpDir == temp1Dir[j])
        ImpComp[k,6] = Imputed[i,"Direction"]  
        k = k+1
      }
    }
  }else{
    for(j in 1:length(temp0Dir)){
      for(l in 1:length(temp1Dir)){
        ImpComp[k,1] = Imputed[i,"Trait_A"]
        ImpComp[k,2] = Imputed[i,"Trait_B"]
        ImpComp[k,3] = (tempImpDir == temp0Dir[j]) == (tempImpDir == temp1Dir[l])
        ImpComp[k,4] = (tempImpDir == temp0Dir[l])
        ImpComp[k,5] = (tempImpDir == temp1Dir[j])
        ImpComp[k,6] = Imputed[i,"Direction"]  
        k = k+1
      }
    }
  }
}

ImpComp = unique(ImpComp)
Issues = ImpComp[which(ImpComp$AllEqual == FALSE),]
Random_F =length(which(Issues$ImputDir == "Random"))
Positive_F = length(which(Issues$ImputDir == "Positive"))
Negative_F = length(which(Issues$ImputDir =="Negative"))

Counts = c(Random_F, Positive_F, Negative_F)
Direction = c("Random", "Positive", "Negative")
Normalized = c((Random_F/nrow(Issues)),(Positive_F/nrow(Issues)), (Negative_F/nrow(Issues)))
Iss_sum = data.frame(Direction, Counts,Normalized)

a = ggplot(data=Iss_sum, aes(x=Direction, y=Normalized)) + geom_bar(stat="identity", fill = "#478A99", width = 0.5)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
ggsave("ImputDiff.pdf")

Match = length(which(ImpComp$AllEqual == TRUE))
Diffs = length(which(ImpComp$AllEqual == FALSE))

DCounts = c(Match, Diffs)
Category = c("Same", "Different")

ImputComp = data.frame(Category, DCounts)
a = ggplot(data=ImputComp, aes(x=Category, y=DCounts)) + geom_bar(stat="identity", fill = "#478A99", width = 0.5)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90))
ggsave("ImputComp.pdf")

rownames(RawData) = RawData[,1]
RawData = RawData[,-1]
RawData=RawData[-which(rownames(RawData) == "Saccharomyces pastorianus"),]
RawData=RawData[-which(rownames(RawData) == "Candida beunavistaensis"),]

rownames(AllImputed) = AllImputed[,1]
AllImputed = AllImputed[,-1]

DataCheck = matrix(nrow = nrow(RawData), ncol = ncol(RawData))
colnames(DataCheck) = colnames(RawData)
rownames(DataCheck) = rownames(RawData)

for(i in 1:nrow(RawData)){
  for(j in 1:ncol(RawData)){
    if(is.na(RawData[i,j]) == F){
      ROW = which(rownames(AllImputed) == rownames(RawData)[i])
      COL = which(colnames(AllImputed) == colnames(RawData)[j])
      DataCheck[i,j] = RawData[i,j] == AllImputed[ROW,COL]
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



