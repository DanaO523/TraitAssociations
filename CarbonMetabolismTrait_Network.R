library("igraph")

CarbonTraits = read.csv("CarbonMetabolismTraits.csv", header = TRUE, check.names = FALSE)
CarbonTraits = CarbonTraits[,1]
AllAsscData = read.csv("PooledAssociation_df.csv", header = TRUE, check.names = F)

Edges = AllAsscData[which(AllAsscData$Trait_A %in% CarbonTraits & AllAsscData$Trait_B %in% CarbonTraits),]

Edges = data.frame(Edges[,c(1,2,10,12)])
TRAITS2 = unique(c(Edges$Trait_A, Edges$Trait_B))

# Edges$EdgeCols = rep("#8A0F31", nrow(Edges))
# Edges[which(Edges$Direction == "Positive"),"EdgeCols"] = "#DAC758"
# Edges = Edges[-which(Edges$Direction == "Shit" |Edges$Direction == "Random"),]

# TRAITS2 = unique(c(Edges$Trait_A, Edges$Trait_B))

# g = graph.data.frame(Edges, directed = FALSE)
# edgeColors = Edges$EdgeCols
# edgeWidth = Edges$AbsDiff/10
# #wc = walktrap.community(g)
# #c = fastgreedy.community(g, weights = edgeWidth)

# c1 = cluster_fast_greedy(g, weights = edgeWidth)
# minC = rep(-Inf, vcount(g))
# maxC = rep(Inf, vcount(g))
# minC[1] <- maxC[1] <- 0

# coords = layout_with_fr(g, minx = minC, maxx = maxC, miny = minC, maxy = maxC)

# layout <- layout.fruchterman.reingold(g)
# pdf(paste("CarbonTrait_AssocGraph_", Sys.Date(), ".pdf", sep = ""), height = 10, width = 10)
# colors <- rainbow(max(membership(fc)))
# plot(c1,g,vertex.color=colors[membership(fc)], vertex.label.color = "black", vertex.label.cex = 1,edge.color=edgeColors, edge.width = edgeWidth,
     # layout=coords, rescale = F, xlim = range(coords[,1]), ylim = range(coords[,2]))

# dev.off()


EdgesPos = Edges[which(Edges$Direction == "Positive"),]

gPos = graph.data.frame(EdgesPos, directed = FALSE)
edgeColors = EdgesPos$EdgeCols
edgeWidth = EdgesPos$AbsDiff/10
wcPos = walktrap.community(gPos)
fcPos = fastgreedy.community(gPos, weights = edgeWidth)

c1 = cluster_fast_greedy(gPos, weights = edgeWidth)
minC = rep(-Inf, vcount(gPos))
maxC = rep(Inf, vcount(gPos))
minC[1] <- maxC[1] <- 0

coords = layout_with_fr(gPos, minx = minC, maxx = maxC, miny = minC, maxy = maxC)

layout <- layout.fruchterman.reingold(gPos)
pdf(paste("CarbonTrait_PosAssocGraph_", Sys.Date(), ".pdf", sep = ""),, height = 10, width = 10, useDingbats = FALSE)
colors <- rainbow(max(membership(fcPos)))
plot(c1,gPos, vertex.label.color = "black", vertex.label.cex = 1,edge.color="black", edge.width = edgeWidth,
     layout=coords, rescale = F, xlim = range(coords[,1]), ylim = range(coords[,2]))

dev.off()


EdgesNeg = Edges[which(Edges$Direction == "Negative"),]

gNeg = graph.data.frame(EdgesNeg, directed = FALSE)
edgeColors = EdgesNeg$EdgeCols
edgeWidth = EdgesNeg$AbsDiff/10
wcNeg = walktrap.community(gNeg)
fcNeg = fastgreedy.community(gNeg)

c1 = cluster_fast_greedy(gNeg, weights = edgeWidth)
minC = rep(-Inf, vcount(gNeg))
maxC = rep(Inf, vcount(gNeg))
minC[1] <- maxC[1] <- 0

coords = layout_with_fr(gNeg, minx = minC, maxx = maxC, miny = minC, maxy = maxC)
#layout <- layout.fruchterman.reingold(gNeg)
pdf(paste("CarbonTrait_NegAssocGraph_", Sys.Date(), ".pdf", sep = ""), height = 10, width = 10, useDingbats = FALSE)
colors <- rainbow(max(membership(fcNeg)))
plot(c1,gNeg, vertex.label.color = "black", vertex.label.cex = 1,edge.color="black", edge.width = edgeWidth,
     layout=coords, rescale = F, xlim = range(coords[,1]), ylim = range(coords[,2]))

dev.off()

Edges_Pos = length(which(Edges$Direction == "Positive"))/nrow(Edges)
Edges_Neg = length(which(Edges$Direction == "Negative"))/nrow(Edges)

Association = c("Positive", "Negative")
Proportion = c(Edges_Pos, Edges_Neg)
CarbAssoc_df = data.frame(Association, Proportion)

a = ggplot(data=CarbAssoc_df, aes(x=Association, y=Proportion, fill = factor(Association))) +geom_bar(stat="identity")+scale_fill_manual(values = c("#8A0F31","#36454F"))
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
ggsave(paste("CarbAssociation_", Sys.Date(), ".pdf", sep = ""))

CarbT_df = data.frame(Trait = character(), Negative = numeric(), Positive = numeric())
for(i in 1:length(TRAITS2)){
  tempDF = Edges[which(Edges$Trait_A == TRAITS2[[i]] | Edges$Trait_B == TRAITS2[[i]]),]
  CarbT_df[i,1] = TRAITS2[[i]]
  CarbT_df[i,2] = length(which(tempDF$Direction == "Negative"))
  CarbT_df[i,3] = length(which(tempDF$Direction == "Positive"))
}
CarbT_df2 = CarbT_df
write.csv(CarbT_df2, file = paste("CarbT_df_", Sys.Date(), ".csv",sep = ""), row.names = F)

CarbT_df2$NegProp = CarbT_df2$Negative/sum(CarbT_df2$Negative)
CarbT_df2$PosProp = CarbT_df2$Positive/sum(CarbT_df2$Positive)
CarbT_df3 = CarbT_df2[,-which(colnames(CarbT_df2) == "Negative" | colnames(CarbT_df2) == "Positive")]

CarbT_df2 = CarbT_df2[order(-CarbT_df2$NegProp, CarbT_df2$PosProp),]
CarbT_df5 = CarbT_df2[,-which(colnames(CarbT_df2) == "NegProp" | colnames(CarbT_df2) ==  "PosProp")]

CarbT_df5 = melt(CarbT_df5)
TROrder = factor(unique(CarbT_df5$Trait))

a = ggplot(data=CarbT_df5, aes(x=Trait, y=value, fill = factor(variable))) +geom_bar(stat="identity", position = "dodge")+scale_fill_manual(values = c("#8A0F31","#36454F"))+scale_x_discrete(limits = TROrder)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
ggsave(paste("CarbAssoc_Count_", Sys.Date(), ".pdf", sep = ""))


TraitCol = "All"
NegProp = CarbAssoc_df[which(CarbAssoc_df$Association == "Negative"),"Proportion"]
PosProp = CarbAssoc_df[which(CarbAssoc_df$Association == "Positive"),"Proportion"]

CarbAssoc_df2 = data.frame(TraitCol, NegProp, PosProp)
colnames(CarbAssoc_df2) = c("Trait", "NegProp", "PosProp")

CarbT_d4 = rbind(CarbT_df3, CarbAssoc_df2)
TraitOrder = CarbT_d4[order(CarbT_d4$NegProp),"Trait"]

CarbT_d4 = melt(CarbT_d4)

#a = ggplot(data=CarbT_d4, aes(x=Trait, y=value, fill = factor(variable))) +geom_bar(stat="identity", position = "dodge")+scale_fill_manual(values = c("#8A0F31","#36454F"))+scale_x_discrete(limits = TraitOrder)
#a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
#ggsave("CarbAssoc_Count.pdf")


FisherTData = AllAsscData[which(AllAsscData$Trait_A %in% CarbonTraits & AllAsscData$Trait_B %in% CarbonTraits),]
Fisher_m = matrix(0,ncol = 2, nrow = 2)
colnames (Fisher_m) = c("Positive", "Negative")
rownames(Fisher_m) = c("Yes", "No")
#FisherTData = FisherTData[-which(FisherTData$Direction == "Random"),]
PosYes = length(which(FisherTData$Direction == "Positive"))
PosNo = length(which(FisherTData$Direction != "Positive"))
NegYes = length(which(FisherTData$Direction == "Negative"))
NegNo = length(which(FisherTData$Direction != "Negative"))

Fisher_m = matrix(c(PosYes,NegYes,PosNo,NegNo), nrow = 2, dimnames = list(Association = c("Positve", "Negative"), Significant = c("Yes", "No")))
Carb_fisher = fisher.test(Fisher_m)
Carb_fisher
save(list = ls(), file = paste("CarbonMetabolismAnalysis_", Sys.Date(), ".rda", sep = ""))



#FisherTest - Enrichments
#Pentoses/Sugar Alcohols Positive
PentOH = matrix(0,nrow =2, ncol = 2)
PentOH[1,1] = 9
PentOH[1,2] = 1
PentOH[2,1] = 2
PentOH[2,2] = 20
PentOH_fisher = fisher.test(PentOH)

# Glucosides Positive
Glucoside = matrix(0,nrow =2, ncol = 2)
Glucoside[1,1] = 7
Glucoside[1,2] = 0
Glucoside[2,1] = 3
Glucoside[2,2] = 22
Glucoside_fisher = fisher.test(Glucoside)

#Contains Galactose
Galactose = matrix(0,nrow =2, ncol = 2)
Galactose[1,1] = 5
Galactose[1,2] = 2
Galactose[2,1] = 0
Galactose[2,2] = 25
Galactose_fisher = fisher.test(Galactose)
