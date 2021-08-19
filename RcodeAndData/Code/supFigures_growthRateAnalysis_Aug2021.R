setwd("C:/Users/mmlam/Desktop/BergmanLabRotation/EcoliExperimentalEvolutionPaper/frontiersInMicrobiology/Code")
# libraries
library(tidyverse)
library(multcomp)
library(ggpubr)
library(latex2exp)
#library(broom)
#library(purrr)
#library(car)
#library(data.table)
library(formattable)
library(cowplot)

# Read the data
ancestral =  read.csv("../Data/growth_ancestrals.csv")
head(ancestral,n=2)

# For each experiment
exp43 = ancestral %>% filter(Temperature == 43)
exp15 = ancestral %>% filter(Temperature == 15)
exp37 = ancestral %>% filter(Temperature == 37)

my_comparisons <- list(c("606P","607P") ,c("REL606","REL607"),c("606P","REL606"),c("607P","REL607"))

## Growth at 43C
# Plot
#pdf("FigureS1_A.pdf")
dev.new()
ggboxplot(exp43, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal=F,
                     label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()

## Growth at 15C
# Plot
#pdf("FigureS1_C.pdf")
dev.new()
ggboxplot(exp15, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",p.adjust.method = "fdr",var.equal
                     =F,label="p.signif",method.args = list(alternative = "greater"))+ylab("Maximum Growth Rate (OD/min)")+ 
  scale_y_continuous(labels = function(x) format(x, scientific = F))+
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()

## Growth at 37C
# Plot
#pdf("FigureS1_B.pdf")
dev.new()
ggboxplot(exp37, x = "Strain", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = my_comparisons,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal
                     =F,label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+ 
  scale_y_continuous(labels = function(x) format(x, scientific = F)) +  
  font("x.text", size = 11)+font("y.text", size = 11)
#dev.off()



#######################################################################
## Figure S2:
# Read the data
dependent =  read.csv("../Data/dependent_15.csv")
dependent$X = factor(paste(dependent$Strain, dependent$Type,sep="_"))
head(dependent,n=2)
table(dependent$X)
levels(dependent$X)
comp_dep = list(c("607P_Acclimation","607P_Condition"), 
                c("606P_Acclimation","606P_Condition"),
                c("F606-2_Acclimation","F606-2_Condition"),
                c("S606-2_Acclimation","S606-2_Condition"),
                c("S607-1_Acclimation","S607-1_Condition"),
                c("S607-2_Acclimation","S607-2_Condition"))
#pdf("FigureS2_dependent_15.pdf")
dev.new()
ggboxplot(dependent, x = "X", y = "Max_deriv",
          color = "Strain", palette = "jco",add="jitter")+
  stat_compare_means(comparisons = comp_dep,method ="t.test",
                     method.args = list(alternative = "greater"),
                     p.adjust.method = "fdr",var.equal
                     =F,label="p.signif")+ylab("Maximum Growth Rate (OD/min)")+ 
  xlab("Strain")+
  scale_y_continuous(labels = function(x) format(x, scientific = F)) +  
  font("x.text", size = 11)+font("y.text", size = 11)+
  theme(axis.text.x=element_text(angle=45, hjust=1))
#dev.off()