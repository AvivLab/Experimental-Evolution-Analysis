#Libraries
library("tidyverse")
library("ggpubr")
library("ComplexHeatmap")

setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")
source("../r_functions/eefe_functions.R")

#Read and prepare the data
data=read.csv("../Data/mutations_populations.csv")
data = data %>% arrange(Population)
data = data %>% mutate(Strain = str_extract(Population,'[0-9][0-9][0-9]'))
head(data)

# Create summary 
resum = as.data.frame(table(data$Population))
resum$Environment  = resum$Var1
resum$Lineage = resum$Var1
#
colnames(resum)[1]="Population"
colnames(resum)[2]="Mutations"
#Recode
resum$Environment = sub(resum$Environment,pattern="F.*",replacement = "Fast")
resum$Environment = sub(resum$Environment,pattern="R.*",replacement = "Random")
resum$Environment = sub(resum$Environment,pattern="S.*",replacement = "Slow")
resum$Lineage = sub(resum$Lineage,pattern="*.606.*",replacement = "606")
resum$Lineage = sub(resum$Lineage,pattern="*.607.*",replacement = "607")
resum$Both = paste(resum$Environment,resum$Lineage,sep=" ")
head(resum)

compare_means(Mutations ~ Both,  data = resum,method = "t.test")

### Pull out evolutionary strategy found via growth analysis and add to mutational data and resum dataframe:
StrategyTable_Updated = read.csv("../Data/August2021_updatedStrategyTable.csv", header = T)
StrategyTable_Updated = rbind(StrategyTable_Updated, c(0, "606P", "Ancestor", 0, 0, "606", "Other"))
StrategyTable_Updated = rbind(StrategyTable_Updated, c(0, "607P", "Ancestor", 0, 0, "607", "Other"))
StrategyTable_Updated$Strategy[grepl("Specialist", StrategyTable_Updated$Strategy, fixed = TRUE)] <- "Specialist"
StrategyTable_Updated$Strategy[grepl("Other", StrategyTable_Updated$Strategy, fixed = TRUE)] <- "NotSignificant"
StrategyTable <- StrategyTable_Updated
data$Strategy = StrategyTable$Strategy[match(data$Population, StrategyTable$Strain)]
resum$Strategy = StrategyTable$Strategy[match(resum$Population, StrategyTable$Strain)]






##############################################################################################
## PDF for number of mutations for each environmental treatment separated by lineage:
#pdf("NumberMutations_boxplot.pdf")
dev.new()
ggboxplot(resum, x = "Both", y = "Mutations",
          color = "Environment", palette = c("#00A08A",'#F98400','#5BBCD6')
          ,add="jitter")+
  ylab("Number of Mutations")+
  xlab("Linage+Environment")+
  font("x.text", size = 11)+font("y.text", size = 11)+theme_xp()
#dev.off()
compare_means(Mutations ~ Both,  data = resum,method = "t.test")


# Same figure/boxplot of number of mutations when normalize each strain to total number of generations:
# import in generations for each strain:
generationsData = read.csv("../Data/populationGenerationsData.csv", header = TRUE)
resum$TotalGenerations = sapply(seq(1:length(resum$Population)),function(x) generationsData$Total.Generations[str_detect(generationsData$ï..Population,toString(resum$Population[x]))])
resum$mutationsNormalizedByTotalGens = resum$Mutations/resum$TotalGenerations
dev.new()
#pdf("NumberOfMutationsByGensForPop_boxplot.pdf")
ggboxplot(resum, x = "Both", y = "mutationsNormalizedByTotalGens",
          color = "Environment", palette = c("#00A08A",'#F98400','#5BBCD6')
          ,add="jitter")+
  ylab("Number of Mutations/Total Gens")+
  xlab("Linage+Environment")+
  font("x.text", size = 11)+font("y.text", size = 11)+theme_xp()
#dev.off()

compare_means(mutationsNormalizedByTotalGens ~ Both,  data = resum,method = "t.test")
t.test(resum$mutationsNormalizedByTotalGens[resum$Environment == "Random"], resum$mutationsNormalizedByTotalGens[resum$Environment == "Fast" | resum$Environment == "Slow"])
t.test(resum$Mutations[resum$Environment == "Random"], resum$Mutations[resum$Environment == "Fast" | resum$Environment == "Slow"])

## PDF for number of mutations for each evolutionary strategy by lineage:
dev.new()
#pdf("NumberMutationsByEvolStrategy_boxplot.pdf")
ggboxplot(resum, x = "Strategy", y = "mutationsNormalizedByTotalGens",
          color = "Environment", palette = c("#00A08A",'#F98400','#5BBCD6')
          ,add="jitter")+
  ylab("Number of Mutations")+
  xlab("Evolutionary Strategy")+
  font("x.text", size = 11)+font("y.text", size = 11)+theme_xp()
#dev.off()



########################################
# Add Environment
data$Environment = data$Population
data$Environment = sub(data$Environment,pattern="F.*",replacement = "Fast")
data$Environment = sub(data$Environment,pattern="R.*",replacement = "Random")
data$Environment = sub(data$Environment,pattern="S.*",replacement = "Slow")
head(data)

# Plot by mutation class
data$Class = factor(data$Class, levels = c("Non_Synonymous", "Small_indel", "Intergenic_snp",
                                           'Stop','Synonymous','Large_deletion',"Large_amplification"))
data$ClassPlot = data$Class

#pdf("BarplotMutationTypesAugust10.pdf")
p <- ggplot(data) + aes(x = Strain, fill = Class)+
  geom_bar(position = "fill") + 
  scale_fill_manual(values= c("#A6CEE3", "#FB9A99",'#FF7F00',
                              '#6A3D9A','#FDBF6F','#E31A1C',"#1F78B4"),
                    labels= c("Nonsynonymous", "Small indel", "Intergenic snv",
                              'Stop','Synonymous','Large deletion',"Large amplification")
  ) + ylab("Relative Abundance") + 
  xlab("Strain")+facet_wrap("Environment")+
  theme(text = element_text(size = 16), legend.position = "bottom", legend.margin = margin(6,5,6,6),
        legend.title = element_blank(),legend.justification=c("left","bottom"))
dev.new()
p
#ggsave("BarplotMutationTypesAugust10.pdf", p, width=20, height=20, 
       #device = "pdf", units = "cm")
#dev.off()


## PDF for different evolution strategies separated by lineage:
#pdf("BarplotMutationTypesByEvolStrategy.pdf")
dev.new()
ggplot(data) + aes(x = Strategy, fill = Class)+
  geom_bar(position = "fill") + 
  scale_fill_manual(values= c("#A6CEE3", "#FB9A99",'#FF7F00',
                              '#6A3D9A','#FDBF6F','#E31A1C',"#1F78B4"),
                    labels= c("Nonsynonymous", "Small indel", "Intergenic snv",
                              'Stop','Synonymous','Large deletion',"Large amplification")
  ) + ylab("Relative Abundance") + 
  xlab("Strategy")+facet_wrap("Strain")+
  theme_xp()
#dev.off()



##############################################
# Plot by mutation categories:
# When consider GroEL/ES as both heath shock and cold shock:
data$Category = factor(data$Category, levels = c("Metabolic","Heat_Shock","Membrane_and_Cell_Wall",
                                                 "Other","Regulatory","Cold_Shock","Transporter","Heat_and_Cold_Shock"),
                       labels = c("Metabolic","Heat Shock","Membrane/Cell Wall",
                                  "Other","Regulatory","Cold Shock","Transporter","Heat/Cold Shock"))
data$CategoryPlot = data$Category
#pdf("BarplotMutationCategoriesForPopsWithHeat_cold_August10.pdf")
p <- ggplot(data) + aes(x = Strain, fill = Category) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values= c("Metabolic" = "cadetblue2", "Heat Shock" = "brown1","Membrane/Cell Wall" = "coral",
                              "Other" = "blueviolet","Regulatory" = "gold","Cold Shock" = "dodgerblue",
                              "Transporter" = "springgreen","Heat/Cold Shock" = "violetred2")#,
                    #labels= c("Metabolic","Heat Shock","Membrane and Cell Wall","Other","Regulatory","Cold Shock","Transporter","Both Shock")
  ) + ylab("Relative Abundance") + 
  xlab("Strain")+facet_wrap("Environment")+
  theme(text = element_text(size = 16), legend.position = "bottom", legend.margin = margin(6,5,6,6),
        legend.title = element_blank(),legend.justification=c("left","bottom"))
dev.new()
p
#ggsave("BarplotMutationCategoriesForPopsWithHeat_cold_August10.pdf", p, width=20, height=20, 
       #device = "pdf", units = "cm")
#dev.off()
########################################################################



###########################################################################################################
## Mann Whitney U one sided test for mutation types (class) like non-synonymous
# Perform below for both "less" and "greater"
# when tied error, build table manually (tied errors are not significant)
# Different mutation classes p-values:
mutCatTable = setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("mutationType", "LineagePval", "EnvPval", "GeneralPval", "SpecialPval","GenVsSpecPval","NotSigPval"))
mutationCatData = as.data.frame(table(data$Class, data$Population))
hypothesisTest = "greater" #less greater or two.sided
for (i in 1:length(unique(data$Class))) {
  mutationCatVar = as.character(unique(data$Class)[i])
  mutTypeOnly = as.data.frame(t(mutationCatData[mutationCatData$Var1==mutationCatVar,3]))
  colnames(mutTypeOnly) = mutationCatData[mutationCatData$Var1==mutationCatVar,2]
  # For Random to Deterministic
  m1<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("R"))), as.numeric(dplyr::select(mutTypeOnly,!contains("R"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For 606 vs 607
  m2<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("606"))), as.numeric(dplyr::select(mutTypeOnly,contains("607"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For generalist vs all specialist and other
  colnames(mutTypeOnly) <- paste(colnames(mutTypeOnly),StrategyTable$Strategy[match(colnames(mutTypeOnly),StrategyTable$Strain)])
  m3<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Generalist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Specialist vs generalist and other
  m4<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Specialist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Generalist vs specialist
  m5<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Not Significant vs generalist & specialists
  m6<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("NotSignificant"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("NotSignificant"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  mutCatTable[i,] = c(mutationCatVar, m2$p.value, m1$p.value, m3$p.value, m4$p.value, m5$p.value, m6$p.value)
}
# 606 lineage has more large deletions than 607: p-value = 0.01898
# Random strains have more Stop mutations that deterministic strains: p-value = 0.0924
# Random strains have less small_indel mutations than deterministic strains: p-value = 0.0505

#######################################################################################
## Diferent Mutation Categories (like heat shock, cold shock, etc): Figure 3D:
# perform below for both "less" and "greater"
# when tied error, build table manually (tied errors are not significant)
# Different mutation classes p-values:
mutCatTable = setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("mutationCategory", "LineagePval", "EnvPval", "GeneralPval", "SpecialPval","GenVsSpecPval","NotSigPval"))
mutationCatData = as.data.frame(table(data$Category, data$Population))
hypothesisTest = "less" #less greater or two.sided
for (i in 1:length(unique(data$Category))) {
  mutationCatVar = as.character(unique(data$Category)[i])
  mutTypeOnly = as.data.frame(t(mutationCatData[mutationCatData$Var1==mutationCatVar,3]))
  colnames(mutTypeOnly) = mutationCatData[mutationCatData$Var1==mutationCatVar,2]
  # For Random to Deterministic
  m1<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("R"))), as.numeric(dplyr::select(mutTypeOnly,!contains("R"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For 606 vs 607
  m2<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("606"))), as.numeric(dplyr::select(mutTypeOnly,contains("607"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # For generalist vs all specialist and other
  colnames(mutTypeOnly) <- paste(colnames(mutTypeOnly),StrategyTable$Strategy[match(colnames(mutTypeOnly),StrategyTable$Strain)])
  m3<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Generalist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Specialist vs generalist and other
  m4<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("Specialist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Generalist vs specialist
  m5<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("Generalist"))), 
                  as.numeric(dplyr::select(mutTypeOnly,contains("Specialist"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  # Not Significant vs generalist & specialists
  m6<-wilcox.test(as.numeric(dplyr::select(mutTypeOnly,contains("NotSignificant"))), 
                  as.numeric(dplyr::select(mutTypeOnly,!contains("NotSignificant"))), paired=FALSE, exact=FALSE, conf.int=TRUE, alternative = hypothesisTest)
  mutCatTable[i,] = c(mutationCatVar, m2$p.value, m1$p.value, m3$p.value, m4$p.value, m5$p.value, m6$p.value)
}

# 606 lineage strains have less regulatory mutations than 607 strains: p-value = 0.0213
# 606 lineage strains have less cold shock mutations than 607 strains: p-value = 0.0979
# 606 lineage strains have more heat shock mutations than 607 strains: p-value = 0.00194
# 606 lineage strains have more transporter mutations than 607 strains: p-value = 0.07082
# Random strains have less regulatory mutations that deterministic strains: p-value = 0.094
# Random strains have less heat and cold shock mutations than deterministic strains: p-value = 0.07089
# Random strains have more membrane and cell wall mutations than deterministic strains: p-value = 0.04522


## Generalist vs specialists
strategyAllData = t(table(resum$Strategy, resum$Population))
fisher.test(cbind(c(0,8),c(9,7)))


###########################################
### Oncoprint:
## Prepare the input
prepared = data %>% group_by(.dots=c("Population","Gene")) %>% mutate(collapsed=paste(Class, collapse = ';'))
prepared = prepared %>% dplyr::select(Gene, Population, collapsed)
prepared = prepared %>% distinct()
prepmat = prepared %>% spread(Population,collapsed,fill="",drop=FALSE)
## Prepare the matrix
mat=as.matrix(prepmat)
rownames(mat) = mat[,1]
mat = mat[, -1]
#Order the matrix
logimat=!(mat == "")
#mat=mat[order(rowSums(logimat),decreasing=T),]
#Save the mutation
prepared$state="1"
mut_mat = prepared %>% spread(Population,state,fill="0",drop=T)
head(mut_mat)

#Get the populations information
pheno = read.csv("../Data/sample_annotations.csv")
Environment=pheno$Environment
Linage = pheno$Linage
ha = HeatmapAnnotation(Environment = Environment, Linage= Linage,
                       col = list(Linage = c("606"='#446455',"607"='#C7B19C'),
                                  Environment = c("Fast" = "#00A08A", "Random" = '#F98400',
                                                  "Slow"='#5BBCD6')
                       ),
                       annotation_height = unit(c(5, 5, 15), "mm"),
                       annotation_legend_param = list(legend_position = "bottom",Environment = list(title = "Environment"),Linage = list(title = "Linage")))

# Orders
StrainOrder = c(grep(colnames(mat), pattern = "R607"),
                grep(colnames(mat), pattern = "F607"),
                grep(colnames(mat), pattern = "S607"),
                grep(colnames(mat), pattern = "R606"),
                grep(colnames(mat), pattern = "F606"),
                grep(colnames(mat), pattern = "S606"))

#OncoPrint
#pdf("oncoPrintFigure3C.pdf")
dev.new()
oncoPrint(mat, column_order = StrainOrder,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE,
          row_names_gp = gpar(fontsize = 7),
          heatmap_legend_param = list(title = "Mutation Classes", 
                                      at = c("Intergenic_snp",'Synonymous',"Non_Synonymous", 'Stop', "Small_indel" ,
                                             'Large_deletion',"Large_amplification") , 
                                      labels = c("Intergenic SNP",'Synonymous',
                                                 "Nonsynonymous",'Stop', 
                                                 "Small Indels",'Large Deletion',
                                                 'Large Amplification')),
          bottom_annotation = ha, pct_gp = gpar(fontsize = 7),
          split = sample(c("606","607"), nrow(mat), replace = TRUE)
          
)
#dev.off()