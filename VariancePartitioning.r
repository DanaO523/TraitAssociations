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
library("Hmisc")
library("vegan")

IsoTraitComb_m = read.csv("IsoTrait_m_2017-01-30.csv", header = TRUE, check.names = F, row.names = 1)
Assoc_df = read.csv("PooledAssociation_df.csv", header = TRUE)
IsoTAssoc_df = read.csv("Assoc_IxT_df.csv", header = TRUE)
BioProp = read.csv("BiologicalProperties.csv", header = TRUE, row.names = 1)
traits = IsoTraitComb_m[,1:48]
load("PGLS.rda")

Isolation_m = IsoTraitComb_m[,49:ncol(IsoTraitComb_m)]
ISOLATIONS = colnames(IsoTraitComb_m[,49:ncol(IsoTraitComb_m)])

##### Isolation Variance Partitioning Based on Traits ####
Iso_mT = t(Isolation_m)

TRAITS = colnames(traits)[c(8:17,19:41,43)]
IsoCarb_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoCarb_m) = TRAITS
rownames(IsoCarb_m) = ISOLATIONS

IsoCarbProp_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoCarbProp_m) = TRAITS
rownames(IsoCarbProp_m) = ISOLATIONS

for(iso in 1:length(ISOLATIONS)){
	TEMP = IsoTraitComb_m[which( IsoTraitComb_m[,which(colnames(IsoTraitComb_m) == ISOLATIONS[[iso]])] == 1),1:48]
	for(tr in 1:length(TRAITS)){
		COL=which(colnames(TEMP)  == TRAITS[tr])
		IsoCarb_m[iso,tr]= length(which(TEMP[,COL] == 1))
		IsoCarbProp_m[iso,tr]= length(which(TEMP[,COL] == 1))/nrow(TEMP)
	}
}

TRAITS = colnames(traits)[45:ncol(traits)]
IsoTemp_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoTemp_m) = TRAITS
rownames(IsoTemp_m) = ISOLATIONS

IsoTempProp_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoTempProp_m) = TRAITS
rownames(IsoTempProp_m) = ISOLATIONS

for(iso in 1:length(ISOLATIONS)){
	TEMP = IsoTraitComb_m[which( IsoTraitComb_m[,which(colnames(IsoTraitComb_m) == ISOLATIONS[[iso]])] == 1),1:48]
	for(tr in 1:length(TRAITS)){
		COL=which(colnames(TEMP)  == TRAITS[tr])
		IsoTemp_m[iso,tr]= length(which(TEMP[,COL] == 1))
		IsoTempProp_m[iso,tr]= length(which(TEMP[,COL] == 1))/nrow(TEMP)
	}
}


TRAITS = colnames(traits)[1:7]
IsoFerm_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoFerm_m) = TRAITS
rownames(IsoFerm_m) = ISOLATIONS

IsoFermProp_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoFermProp_m) = TRAITS
rownames(IsoFermProp_m) = ISOLATIONS

for(iso in 1:length(ISOLATIONS)){
	TEMP = IsoTraitComb_m[which( IsoTraitComb_m[,which(colnames(IsoTraitComb_m) == ISOLATIONS[[iso]])] == 1),1:48]
	for(tr in 1:length(TRAITS)){
		COL=which(colnames(TEMP)  == TRAITS[tr])
		IsoFerm_m[iso,tr]= length(which(TEMP[,COL] == 1))
		IsoFermProp_m[iso,tr]= length(which(TEMP[,COL] == 1))/nrow(TEMP)
	}
}


TRAITS = colnames(traits)[c(18,42,44)]
IsoMisc_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoMisc_m) = TRAITS
rownames(IsoMisc_m) = ISOLATIONS

IsoMiscProp_m = matrix(0, nrow = length(ISOLATIONS), ncol = length(TRAITS))
colnames(IsoMiscProp_m) = TRAITS
rownames(IsoMiscProp_m) = ISOLATIONS

for(iso in 1:length(ISOLATIONS)){
	TEMP = IsoTraitComb_m[which( IsoTraitComb_m[,which(colnames(IsoTraitComb_m) == ISOLATIONS[[iso]])] == 1),1:48]
	for(tr in 1:length(TRAITS)){
		COL=which(colnames(TEMP)  == TRAITS[tr])
		IsoMisc_m[iso,tr]= length(which(TEMP[,COL] == 1))
		IsoMiscProp_m[iso,tr]= length(which(TEMP[,COL] == 1))/nrow(TEMP)
	}
}

IsoVarPart = varpart(Iso_mT, IsoCarb_m,IsoTemp_m, IsoFerm_m,IsoMisc_m)

pdf("IsolationVarianceParition.pdf")
plot(IsoVarPart, cutoff = 0.0005, digits = 2)
dev.off()


#### Trait Association Variance ####
TRAITS = colnames(traits)

TRCOMB = list()
k = 1
for(i in 1:(length(TRAITS)-1)){
	for(j in (i+1):length(TRAITS)){
		TRCOMB[k] = paste(TRAITS[[i]], TRAITS[[j]], sep = "|")
		k = k+1
	}
}
TRCOMB = unlist(TRCOMB)


PresCombTrait_m = matrix(nrow = nrow(traits), ncol = length(TRCOMB))
rownames(PresCombTrait_m)  = rownames(traits)
colnames(PresCombTrait_m) =  TRCOMB

for(i in 1:(length(TRAITS)-1)){
	for(j in (i+1):length(TRAITS)){
		COLS=which(colnames(traits) == TRAITS[[i]] | colnames(traits) == TRAITS[[j]])
		TEMP = traits[,COLS]
		TEMPSUMS = matrix(rowSums(TEMP))
		TEMPSUMS[which(TEMPSUMS == 1),] = 0
		TEMPSUMS[which(TEMPSUMS == 2),] = 1
		TRCOLS = which(colnames(PresCombTrait_m) == paste(TRAITS[[i]], TRAITS[[j]], sep = "|"))
		PresCombTrait_m[,TRCOLS] = TEMPSUMS
	}
}

BIOPROPS = colnames(BioProp)
BioPropTrait_m = matrix(nrow = nrow(traits), ncol = length(BIOPROPS))
rownames(BioPropTrait_m)  = rownames(traits)
colnames(BioPropTrait_m) =  BIOPROPS

PropBioPropTrait_m = matrix(nrow = nrow(traits), ncol = length(BIOPROPS))
rownames(PropBioPropTrait_m)  = rownames(traits)
colnames(PropBioPropTrait_m) =  BIOPROPS


for(i in 1:nrow(traits)){
	for(j in 1:length(BIOPROPS)){
		TEMP = rownames(BioProp)[which(BioProp[,j] == 1)]
		DATA = data.frame(traits[i,which(colnames(traits) %in% TEMP)], check.names = F)
		BioPropTrait_m[i,j] = length(which(DATA[1,] == 1))
		PropBioPropTrait_m[i,j] = length(which(DATA[1,] == 1))/length(TEMP)
	}
}

Phy_m = cophenetic(rootedtree)

COLUMNS = colnames(Phy_m)
COL_ed = lapply(COLUMNS, function(y) gsub("_", " ", y)) 
COL_ed = unlist(COL_ed)

colnames(Phy_m) = COL_ed

ROWS = rownames(Phy_m)
ROW_ed = lapply(ROWS, function(y) gsub("_", " ", y))
ROW_ed = unlist(ROW_ed)

rownames(Phy_m) = ROW_ed

ROWKEEP =which(rownames(PresCombTrait_m) %in% rownames(Phy_m))
PresCom_limited = PresCombTrait_m[ROWKEEP,]


ROWORDER = rownames(PresCom_limited)
OrderedPhy_m = matrix(0, nrow = length(ROWORDER), ncol = ncol(Phy_m))

for(i in 1:length(ROWORDER)){
OrderedPhy_m[i,] =  Phy_m[which(rownames(Phy_m) == ROWORDER[[i]]),]
}
rownames(OrderedPhy_m) = ROWORDER
colnames(OrderedPhy_m) = colnames(Phy_m)

OrderedPhy_m2=OrderedPhy_m[which(colnames(OrderedPhy_m) %in% ROWORDER)]

Iso_limited  = Isolation_m[which(rownames(Isolation_m) %in% rownames(PresCom_limited)),]
Bio_limited = BioPropTrait_m[which(rownames(BioPropTrait_m) %in% rownames(PresCom_limited)),]

VarPart_limited = varpart(PresCom_limited, Iso_limited,OrderedPhy_m2,Bio_limited)
pdf("VarianceParition_limited.pdf")
plot(VarPart_limited, cutoff = 0.000000005, digits = 5)
dev.off()

test = rda(PresCom_limited, Iso_limited,OrderedPhy_m2,Bio_limited)
