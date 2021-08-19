#Libraries
library('tidyverse')
library("rlang")
library("dplyr")

######################################################
# Read the results:
setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Data")
files = list.files(pattern = "*BothDirectionality.csv")
setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")
direccio = read.csv(paste("../Data/",files[1],sep=""))

for (i in 2:length(files)){
  directe = read.csv(paste("../Data/",files[i],sep="")) 
  direccio = rbind(direccio,directe)
}
head(direccio)

# Remove Dim 4 from analysis:
removedPc4 <- direccio[direccio$PC != "Dim.4",]

#Update strategy for what found on August 10th with Kenny Ye
StrategyTable_Updated = read.csv("../Data/August2021_updatedStrategyTable.csv", header = T)
StrategyTable_Updated = rbind(StrategyTable_Updated, c(0, "606P", "Ancestor", 0, 0, "606", "Other"))
StrategyTable_Updated = rbind(StrategyTable_Updated, c(0, "607P", "Ancestor", 0, 0, "607", "Other"))
StrategyTable_Updated$Strategy[grepl("Specialist", StrategyTable_Updated$Strategy, fixed = TRUE)] <- "Specialist"
StrategyTable_Updated$Strategy[grepl("Other", StrategyTable_Updated$Strategy, fixed = TRUE)] <- "Not Significant"

removedPc4$Strategy <- StrategyTable_Updated$Strategy[match(removedPc4$Name, StrategyTable_Updated$Strain)]

################################################3
## Condense to 5 categories for hug & gaut analysis: 
# 1. restorative (over-restored, restored, and partially restored), 
# 2. unrestorative (unrestored), 
# 3. novel, 
# 4. reinforced, 
# 5. other (uninformative and NA).

## Function to condense categories of restoration analysis results:
direccio <- removedPc4
direccio$CondensedComparison = direccio$comparison
direccio$CondensedComparison = sub(direccio$CondensedComparison,pattern="Partially Restored",replacement = "Restorative")
direccio$CondensedComparison = sub(direccio$CondensedComparison,pattern="Over-restored",replacement = "Restorative")
direccio$CondensedComparison = sub(direccio$CondensedComparison,pattern="Restored",replacement = "Restorative")
direccio$CondensedComparison = sub(direccio$CondensedComparison,pattern="Uninformative",replacement = "Uninformative")


#pdf("Figure4D_updatedAugust2021.pdf")
plotDataFrame <- direccio[direccio$Condition==15,] 
dev.new()
ggplot(plotDataFrame) + aes(x = PC, fill = CondensedComparison)+
  geom_bar(position="fill") +
  ylab("Relative Abundance") + 
  xlab("Environment")+facet_grid(Treatment~referencia~Strategy)+ theme(text = element_text(size = 14),axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()

# For each biolog temperature treated, 15 and 43C, by strain (606 vs 607) and treatment (fast, slow, or random):
dev.new()
ggplot(direccio) + aes(x = PC, fill = CondensedComparison)+
  geom_bar(position = "fill") +
  ylab("Relative Abundance") + 
  xlab("Environment")+facet_grid(Strain~referencia~Treatment)

# Do fishers exact test to look for any significance in signature of 15C strain vs 43C strains: 
restorationData15 = as.data.frame(table(direccio$CondensedComparison[direccio$Condition == 15], direccio$referencia[direccio$Condition == 15]))
restorationData43 = as.data.frame(table(direccio$CondensedComparison[direccio$Condition == 43], direccio$referencia[direccio$Condition == 43]))
responseType = "Unrestored" # "Unrestored"    "Restorative"   "Reinforced"    "Uninformative" "Novel"
fisher.test(cbind(
  c(restorationData15$Freq[restorationData15$Var1 == responseType & restorationData15$Var2 == "15"],
    sum(restorationData15$Freq[restorationData15$Var1 != responseType & restorationData15$Var2 == "15"])),
  c(restorationData43$Freq[restorationData43$Var1 == responseType & restorationData43$Var2 == "43"],
    sum(restorationData43$Freq[restorationData43$Var1 != responseType & restorationData43$Var2 == "43"]))))


# Do one-sided Mann Whitney U test for both temperatures to look for any significance between lineage, deterministic vs random treatment (treatment), and generalist vs specialist:
restorationTable = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("restorationCategory", "LineagePval", "EnvPval", "GeneralPval", "SpecialPval","GenVsSpecPval"))
restorationData = as.data.frame(table(direccio$CondensedComparison[(direccio$referencia == 15 | direccio$referencia == 43) 
                                                                   & direccio$Condition == 15], 
                                      direccio$Name[(direccio$referencia == 15 | direccio$referencia == 43) 
                                                    & direccio$Condition == 15]))
hypothesisTest = "greater" # do for both "less" and "greater"
for (i in 1:length(unique(restorationData$Var1))) {
  mutationCatVar = as.character(unique(restorationData$Var1)[i])
  mutTypeOnly = as.data.frame(t(restorationData[restorationData$Var1==mutationCatVar,3]))
  colnames(mutTypeOnly) = restorationData[restorationData$Var1==mutationCatVar,2]
  # For Random to Deterministic
  m1<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("R"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("R"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For 606 vs 607
  m2<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("606"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("607"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For generalist vs all specialist and other
  colnames(mutTypeOnly) <- paste(colnames(mutTypeOnly),StrategyTable_Updated$Strategy[match(colnames(mutTypeOnly),StrategyTable_Updated$Strain)])
  m3<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Generalist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Specialist vs generalist and other
  m4<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Specialist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Generalist vs specialist
  m5<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  restorationTable[i,] = c(mutationCatVar, m2$p.value, m1$p.value, m3$p.value, m4$p.value, m5$p.value)
}


# Do one-sided Mann Whitney U test for 15C OR 43C only to look for any significance between lineage, deterministic vs random treatment (treatment), and generalist vs specialist:
tempLookingAt = 43 # 15 or 43
restorationTable = setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("restorationCategory", "LineagePval", "EnvPval", "GeneralPval", "SpecialPval","GenVsSpecPval"))
restorationData = as.data.frame(table(direccio$CondensedComparison[direccio$referencia == tempLookingAt & direccio$Condition == tempLookingAt], 
                                      direccio$Name[direccio$referencia == tempLookingAt & direccio$Condition == tempLookingAt]))
hypothesisTest = "less" # do for both "less" and "greater"
for (i in 1:length(unique(restorationData$Var1))) {
  mutationCatVar = as.character(unique(restorationData$Var1)[i])
  mutTypeOnly = as.data.frame(t(restorationData[restorationData$Var1==mutationCatVar,3]))
  colnames(mutTypeOnly) = restorationData[restorationData$Var1==mutationCatVar,2]
  # For Random to Deterministic
  m1<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("R"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("R"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For 606 vs 607
  m2<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("606"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("607"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For generalist vs all specialist and other
  colnames(mutTypeOnly) <- paste(colnames(mutTypeOnly),StrategyTable_Updated$Strategy[match(colnames(mutTypeOnly),StrategyTable_Updated$Strain)])
  m3<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Generalist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Specialist vs generalist and other
  m4<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Specialist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Generalist vs specialist
  m5<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  restorationTable[i,] = c(mutationCatVar, m2$p.value, m1$p.value, m3$p.value, m4$p.value, m5$p.value)
}
