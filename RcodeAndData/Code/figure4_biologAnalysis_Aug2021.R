#Libraries
library("tidyverse")
library("ggpubr")
library("pheatmap")
library("ComplexHeatmap")
library("FactoMineR")
library("factoextra")
library("UpSetR")
library("pca3d")
library("gplots")
library("cowplot")
library("latex2exp")
library("ggcorrplot")
library("magrittr")

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")
source("../r_functions/eefe_functions.R")

## PCA of biolog data:
#Read and prepare the data
data606 = read.csv("../Data/normalized606Wells86_94RemovedBiologData.csv")
data607 = read.csv("../Data/normalized607Wells86_94RemovedBiologData.csv")
wells = read.csv("../Data/gen_biolog.csv")
wells =  wells %>% filter(ï..Well %in% colnames(data606))
# Read the strategies file
strategy = read.csv("../Data/StrategyTable.csv") %>% dplyr::select(Name, Strategy)
data606 = inner_join(data606,strategy, by = "Name") %>% dplyr::select(1:5,Strategy,everything())
data607 = inner_join(data607,strategy, by = "Name") %>% dplyr::select(1:5,Strategy,everything())
colnames(data606)[7:98] = as.character(wells$Assay)
colnames(data607)[7:98] = as.character(wells$Assay)


#################################################################################################
#################################################################################################
## Both 606P and 607P strains together:
allData <- rbind(data606,data607)
# If want ancestors separated by lineage:
#allData$Merged =  paste(allData$Strain,allData$Treatment,allData$Condition,sep="") 
# If do NOT want ancestors separated by lineage, so only by treatment and condition:
allData$Merged =  paste(allData$Treatment,allData$Condition,sep="") 
allData = allData %>% dplyr::select(Merged,everything())
ColFactor =  mutate(allData,Colors = ifelse(Treatment == "Ancestor", Merged, Condition))
# Log transform
dat = allData[,-c(1:7)]
colnames(dat) = wells$Assay
rownames(dat) = make.unique(paste(allData$Name,allData$Condition,sep="-"),sep="-")
logdat = log1p(dat)
# PCA
res.pca = prcomp(logdat, scale = F)
# For writing scores from PCA:
res.pca.pca = PCA(logdat, scale.unit = F, graph = F)
scores =cbind(allData[,c(1:7)],res.pca.pca$ind$coord)
PCsToPlot = c(1,2)
# PCA plot with 606 vs 607 labeled as different shapes and then colored by treatment
# ****Ellipses should be for combined 606 and 607 ancestors****
# If want ancestors separated by lineage:
#colorByBiologCondition = list(TreatmentTemp = c("15"="navy","43" = "firebrick", "37" = 'orange',"606Ancestor15"="grey85","606Ancestor43"="grey86","606Ancestor37"="grey87","607Ancestor15"="grey52","607Ancestor37"="grey87","607Ancestor43"="grey53"))
# If want ancestors NOT separated by lineage:
colorByBiologCondition = list(TreatmentTemp = c("15"="navy","43" = "firebrick", "37" = 'orange',"Ancestor15"="grey53","Ancestor43"="grey52","Ancestor37"="grey51"))
#pdf("Figure4B_updatedAugust2021.pdf")
p <- fviz_pca_ind(res.pca,
             axes = PCsToPlot,
             geom.ind = "point", # show points only (but not "text")
             col.ind = "white", # color by groups
             palette = c("navy","orange","firebrick","grey53","grey52","grey51"), 
             addEllipses = T, 
             ellipse.type = "norm",
             label = "var",
             col.var = "white", repel = TRUE,
             legend.title = "Factor",
             pointshape = 21,
             pointsize = 1.5,
             fill.ind = as.character(ColFactor$Colors),  
             ellipse.level = 0.8
)
p <- fviz_add(p, res.pca$x[(allData$Strain=="607"),], axes = PCsToPlot, 
         color = colorByBiologCondition$TreatmentTemp[ColFactor$Colors][allData$Strain=="607"],
         shape = 8, addlabel = FALSE) 
fviz_add(p, res.pca$x[(allData$Strain=="606"),], axes = PCsToPlot, 
         color = colorByBiologCondition$TreatmentTemp[ColFactor$Colors][allData$Strain=="606"],
         shape = 15, addlabel = FALSE)
#dev.off()
dev.new()
p


# Save the scores:
#write.csv(scores,"../Data/both606And607_scores.csv")
# Talus plot for the informative principal components for Figure s10:
eigenvalues = res.pca.pca$eig[,1]
cut = 3
########## RUN THIS HERE ###########
#pdf("Talus_Both606and607.pdf")
TalPlot(eigenvalues,cut)
#dev.off()
####################################



#################################
#Figure s10:
p1 = fviz_contrib(res.pca.pca, choice = "var", axes =1, top = 10) + theme(axis.text=element_text(size=6))
p2 = fviz_contrib(res.pca.pca, choice = "var", axes =2, top = 10) + theme(axis.text=element_text(size=6))
p3 = fviz_contrib(res.pca.pca, choice = "var", axes =3, top = 10) + theme(axis.text=element_text(size=6))
#p4 = fviz_contrib(res.pca, choice = "var", axes =4, top = 10) + theme(axis.text=element_text(size=6))
#pdf("Both606and607_Contributions.pdf")
plot_grid(p1, p2, p3,labels = c('', '',""))
#dev.off()
vint1 = p1$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint2 = p2$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vint3 = p3$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
#vint4 = p4$data %>% arrange(-contrib) %>% top_n(contrib,n=10)
vintmes = unique(c(as.character(vint1$name),as.character(vint2$name),
                   as.character(vint3$name)))
correlacions = res.pca.pca$var$cor[vintmes,1:3]
# Plot
#pdf("Both606and607_TopCorrelation3PCs.pdf")
ggcorrplot(t(correlacions), tl.cex  = 10)
#dev.off()


